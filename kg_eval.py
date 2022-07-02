import re
from collections import defaultdict, Counter
import pandas as pd


def getGeneName(s):
    return s.split("*")[0]


def extractID(name):
    return re.findall(r"\.(\w+)\.read", name)[0]


class EvaluateKIR:
    def __init__(self, summary_csv: str):
        """ Get answer alleles from summary.csv """
        data = pd.read_csv(summary_csv, sep="\t", dtype=str)
        print(data)
        # this is ordered dict
        self.ans = {i.id: sorted(i.alleles.split("_")) for i in data.itertuples()}

    def getAns(self, id):
        """ Return allelelist of str """
        return self.ans[id]

    def getAnsCN(self, id):
        alleles = self.getAns(id)
        return Counter(map(getGeneName, alleles))

    def compareSample(self, id, predict_list):
        return EvaluateKIR.test(self.getAns(id), predict_list, title=f"Sample {id}")

    def compareCohert(self, predict_cohert, skip_empty=False):
        """
        Compare called alleles from all samples

        Args:
          predict_cohert(Any): It can be dict[str, list[str]] or list[list[str]]

        Return:
          summary(dict[str, int]): The number of each metrics
        """
        # if list -> compare the ans with same order
        if type(predict_cohert) is list:
            ids = list(self.ans.keys())
            if skip_empty:
                ids = ids[:min(len(ids), len(predict_cohert))]
            results = [self.compareSample(id, predict_cohert[i]) for i, id in enumerate(ids)]
        # if list -> compare the ans with same id
        elif type(predict_cohert) is dict:
            ids = self.ans.keys()
            if skip_empty:
                ids =  ids & predict_cohert.keys()
            ids = sorted(ids)
            results = [self.compareSample(id, predict_cohert.get(id, [])) for id in ids]
        else:
            raise NotImplementedError
        return EvaluateKIR.summarize(results)

    @staticmethod
    def test(ans_list, predict_list, title=""):
        """ Compare two alleles set """
        predict_list = sorted(predict_list)
        comparison_tuple = []

        for gene in set(getGeneName(i) for i in ans_list + predict_list):
            # if gene in ["KIR2DP1", "KIR3DP1"]:
            #     continue
            if "KIR2DL5*unresolved" in predict_list and gene in ["KIR2DL5A", "KIR2DL5B"]:
                continue  # skip 2DL5A 2DL5B only run when gene == "KIR2DL5"
            comparison_tuple.extend(
                EvaluateKIR.testPerGene([i for i in ans_list if i.startswith(gene)],
                                        [i for i in predict_list if i.startswith(gene)])
            )
        comparison_tuple.sort(key=lambda i: i[1] or i[2])
        print(title)
        for t, a, b in comparison_tuple:
            if t == "Match7":
                print(f"{a:18} OK {b:18}")
            elif t == "Match5":
                print(f"{a:18} <5 {b:18}")
            elif t == "Match3":
                print(f"{a:18} <3 {b:18}")
            elif t == "MatchGene":
                print(f"{a:18} <0 {b:18}")
            elif t == "FN":
                print(f"{a:18} <-")
            elif t == "FP":
                print(f"{'':18} -> {b:18}")
        return comparison_tuple

    @staticmethod
    def testPerGene(a_list, b_list):
        """
        Compare two alleles set

        (All alleles in these two set must in the same gene)

        Args:
          a_list: list[str]
          b_list: list[str]

        Return:
          list[ tuple[str, str, str] ]: list of comparison
            each comparison contains three items
            * type
            * a_allele
            * b_allele
        """
        # Find perfect match
        for allele in list(b_list):
            if allele in a_list:
                a_list.remove(allele)
                b_list.remove(allele)
                yield ("Match7", allele, allele)

        # Find match 5 digits
        for allele_b in list(b_list):
            for allele_a in a_list:
                if allele_a[:allele_a.index("*") + 6] == allele_b[:allele_b.index("*") + 6]:
                    a_list.remove(allele_a)
                    b_list.remove(allele_b)
                    yield ("Match5", allele_a, allele_b)
                    break

        # Find match 3 digits
        for allele_b in list(b_list):
            for allele_a in a_list:
                if allele_a[:allele_a.index("*") + 4] == allele_b[:allele_b.index("*") + 4]:
                    a_list.remove(allele_a)
                    b_list.remove(allele_b)
                    yield ("Match3", allele_a, allele_b)
                    break

        # Find match < 3 digits
        for allele_a, allele_b in zip(list(a_list), list(b_list)):
            yield ("MatchGene", allele_a, allele_b)
            a_list.remove(allele_a)
            b_list.remove(allele_b)

        # copy number error
        for allele in a_list:
            yield ("FN", allele, None)
        for allele in b_list:
            yield ("FP", None, allele)

    @staticmethod
    def summarize(comparison_list):
        """
        Summarize the comparison

        Args:
          comparison_list(list[list[tuple]]): The comparisons of each sample
        Return:
          summary(dict[str, int]): The number of each metrics
        """
        summary = defaultdict(int)
        # TP, FP, FN, match_gene, match_3, match_5, total, cn_error = 0, 0, 0, 0, 0, 0, 0, 0
        for result in comparison_list:
            summary['total'] += len([i[1] for i in result if i[1]])
            summary['total_sample'] += 1
            if all(i[0] == 'FN' for i in result):
                summary['fail_allele'] += len(result)
                summary['fail_sample'] += 1
                continue
            for i in result:
                summary[i[0]] += 1

        print(f"Total alleles       = {summary['total']} (sample = {summary['total_sample']})")
        print(f"  * Cannot_called   = {summary['fail_allele']} (sample = {summary['fail_sample']})")
        print(f"  * TP              = {summary['Match7']}")
        print(f"  * Match_5_digits  = {summary['Match5']}")
        print(f"  * Match_3_digits  = {summary['Match3']}")
        print(f"  * Match_gene      = {summary['MatchGene']}")
        print(f"  * Copy_number_err = {summary['FN'] + summary['FP']}")
        print(f"    * FN = {summary['FN']}")
        print(f"    * FP = {summary['FP']}")
        return summary


if __name__ == "__main__":
    # test
    # result = EvaluateKIR.test(["KIR3DL3*097", "KIR3DL3*0020605", "KIR2DS2*0010112", "KIR2DS2*0010102", "KIR2DL2*0010102", "KIR2DL2*0010101", "KIR2DL5A*00107", "KIR2DL5A*030", "KIR2DL5B*0390102", "KIR2DS3*009", "KIR2DP1*0020105", "KIR2DL1*0030230", "KIR3DP1*026", "KIR3DP1*00312", "KIR2DL4*0010302", "KIR2DL4*045", "KIR3DS1*0130108", "KIR3DS1*0130102", "KIR2DS5*0020105", "KIR2DS5*0020107", "KIR2DS1*013", "KIR2DS1*0020103", "KIR3DL2*0070102", "KIR3DL2*023"], ["KIR3DL3*097", "KIR2DS2*00102", "KIR2DS2*003", "KIR2DL2*0010103"])
    # EvaluateKIR.summarize([result])

    # evaluate PING
    # kir = EvaluateKIR("linnil1_syn_wide/linnil1_syn_wide.summary.csv")

    # 100 samples
    # ping_called_alleles = readPingResult(f"/home/linnil1/kir/PING/PING_test/linnil1_syn_wide_result/finalAlleleCalls.csv")
    # kir.compareCohert(ping_called_alleles)

    # 10 samples
    # ping_called_alleles = {k: v for k, v in ping_called_alleles.items() if k.startswith("0")}
    # kir.compareCohert(ping_called_alleles, skip_empty=True)

    kir = EvaluateKIR("linnil1_syn_30x/linnil1_syn_30x.summary.csv")
    ping_called_alleles = readPingResult(f"./PING/data_linnil1_syn_30x.result/finalAlleleCalls.csv")
    kir.compareCohert(ping_called_alleles)
