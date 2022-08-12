import json
from collections import defaultdict
import pandas as pd

from gk_hisat2 import ReadsAndVariantsData, loadReadsAndVariantsData, PairRead
from gk_multi_allele_typing import AlleleTyping
from gk_hisat2em import Hisat2AlleleResult, printHisatTyping


def removeMultipleMapped(reads_data: ReadsAndVariantsData) -> ReadsAndVariantsData:
    """ actually remove NH != 1 reads """
    return {
        'variants': reads_data['variants'],
        'reads': list(filter(lambda i: i.multiple != 1, reads_data['reads'])),
    }


def groupReads(reads: list[PairRead]) -> dict[str, list[PairRead]]:
    """ Group the reads by gene name (backbone) """
    gene_reads = defaultdict(list)
    for read in reads:
        gene_reads[read.backbone].append(read)
    return gene_reads


class Typing:
    """ Typing allele tools """
    def __init__(self, filename: str):
        raise NotImplementedError

    def typingPerGene(self, gene: str, cn: int) -> list[str]:
        raise NotImplementedError

    def typing(self, gene_cn: dict[str, int]) -> list[str]:
        """ Allele Typing for all genes with CN """
        predict_alleles = []
        for gene, cn in gene_cn.items():
            if not cn:
                continue
            predict_alleles.extend(self.typingPerGene(gene, cn))
        return predict_alleles


class TypingWithPosNegAllele(Typing):
    """ Our proposal method """
    def __init__(self, filename_variant_json: str, multiple=False):
        reads_data = loadReadsAndVariantsData(filename_variant_json)
        if not multiple:
            reads_data = removeMultipleMapped(reads_data)
        self.gene_reads = groupReads(reads_data['reads'])
        self.variants = reads_data['variants']

    def typingPerGene(self, gene: str, cn: int) -> list[str]:
        typ = AlleleTyping(self.gene_reads[gene], self.variants, {})
        res = typ.typing(cn)
        return res.selectBest()


class TypingWithReport(Typing):
    """ Typing alleles from report by abundance value"""
    def __init__(self, filename_report: str):
        self.report_data = json.load(open(filename_report))
        for gene, alleles in self.report_data.items():
            self.report_data[gene] = [Hisat2AlleleResult(**i) for i in alleles]

        printHisatTyping(self.report_data)

    def typingPerGene(self, gene: str, cn: int) -> list[str]:
        """
        Example:
        case1: a1=0.66 a2=0.33
          CN=1: a1
          CN=2: a1 a1
          CN=3: a1 a1 a2
        case2: a1=0.90 a2=0.10
          CN=1: a1
          CN=2: a1 a1
          CN=3: a1 a1 a1
        """
        report = self.report_data[gene]

        est_prob = 1 / cn
        alleles = sorted(report, key=lambda i: -i.prob)
        called_alleles = []
        for allele in alleles:
            pred_count = max(1, round(allele.prob / est_prob))
            for i in range(min(cn, pred_count)):
                called_alleles.append(allele.allele)
            cn -= pred_count
            if cn <= 0:
                break
        return called_alleles


def loadCN(filename_cn: str) -> dict[str, int]:
    """ Load CN data """
    data = pd.read_csv(filename_cn, sep="\t", index_col=[0])
    return data.to_dict()['cn']


if __name__ == "__main__":
    prefix = "data/linnil1_syn_30x.00.index_kir_2100_ab_2dl1s1_muscle_mut01_graph.variant"
    method = "pv"

    if method == "pv":
        filename_read_variant = prefix + ".json"
        t = TypingWithPosNegAllele(filename_read_variant)  # type: Typing
    else:
        filename_em = prefix + ".em.json"
        t = TypingWithReport(filename_em)

    prefix += ".no_multi.depth.cn_p75"
    filename_cn = prefix + ".tsv"
    cn = loadCN(filename_cn)
    print(cn)

    called_alleles = t.typing(cn)
    print(called_alleles)
    df = pd.DataFrame([{
        'name': prefix,
        'alleles': "_".join(called_alleles),
    }])
    filename_result = prefix + ".pn.tsv"
    df.to_csv(filename_result, sep="\t", index=False)
