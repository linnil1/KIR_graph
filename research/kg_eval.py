from enum import Enum
from typing import Iterator, Callable, Iterable, Any
from itertools import chain, repeat
from dataclasses import dataclass
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import re

from Bio import pairwise2, SeqIO, SeqRecord
from Bio.Seq import Seq
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from namepipe import NamePath
from graphkir.plot import showPlot
from graphkir.utils import getGeneName, getAlleleField,  limitAlleleField


CohortAlleles = dict[str, list[str]]  # it's ordered dict too
CohortSequences = dict[str, dict[str, SeqRecord]]  # it's ordered dict too


class MatchType(Enum):
    """
    Possible types of matching two alleles (left = answer, right = predict)

    * MATCH7: answer == predict
    * MATCH5: KIR2DL1*0010101, KIR2DL1*0010103
    * MATCH5: KIR2DL1*0010101, KIR2DL1*00101
    * MATCH3: KIR2DL1*0010101, KIR2DL1*0010201
    * MATCH3: KIR2DL1*00101,   KIR2DL1*00102
    * GENE:   KIR2DL1*001,     KIR2DL1*002
    * FN:     KIR2DL1*0010101, NOT FOUND
    * FP:     NOT FOUND,       KIR2DL1*0010101
    """
    MATCH7    = 0b11111
    MATCH5    = 0b11101
    MATCH3    = 0b11001
    MATCHGENE = 0b10001
    FN        = 0b10000
    FP        = 0b00001
    NONE      = 0b00000


@dataclass
class MatchResult:
    """ The pair of alleles from answer and predict alleles """
    answer_allele: str  # order is important
    predit_allele: str
    answer_allele_full: str  # novel name
    predit_allele_full: str  # novel name
    match_type: MatchType = MatchType.NONE

    # other info
    answer_seq: str = ""
    predit_seq: str = ""
    answer_allele_length: int = 0
    predit_allele_length: int = 0
    base_diff: int = 0

    def __lt__(self, other: 'MatchResult') -> bool:
        return (self.answer_allele or self.predit_allele) \
               < (other.answer_allele or other.predit_allele)


CohortMatchResult = dict[str, list[MatchResult]]  # it's ordered dict too


def groupByGene(alleles: list[str]) -> dict[str, list[str]]:
    """ return dict[gene_name, list[allele_name]] """
    gene_alleles = defaultdict(list)
    for name in alleles:
        gene_alleles[getGeneName(name)].append(name)
    return gene_alleles


def readAnswerAllele(summary_tsv: str) -> CohortAlleles:
    """ Read answer allele """
    data = pd.read_csv(summary_tsv, sep="\t", dtype=str)
    return {i.id: sorted(i.alleles.split("_")) for i in data.itertuples()}


def saveCohortAllele(data: CohortAlleles, summary_tsv: str) -> None:
    """ Read answer allele """
    predict_list = []
    for id, alleles in data.items():
        predict_list.append(
            {
                "id": id,
                "alleles": "_".join(alleles),
                "name": "." + str(id) + ".",
            }
        )
    df = pd.DataFrame(predict_list)
    df.to_csv(summary_tsv, index=False, sep="\t")


def extractID(name: str) -> str:
    """
    Extract the sample id from filename

    Example:
      Input: "linnil1_syn_30x_seed87.00.index.1.02.fa"
      Output: "00"
    """
    return re.findall(r"\.(\d{2})\.", name)[0]


def readPredictResult(tsv_file: str,
                      extract_func: Callable[[str], str] = extractID
                      ) -> CohortAlleles:
    """ Read predict alleles (same format as summary.csv) """
    df = pd.read_csv(tsv_file, sep='\t', dtype=str)
    data = {}
    for i in df.itertuples():
        if "id" in df.columns:
            id = i.id
        else:
            id = extractID(i.name)
        data[id] = sorted(i.alleles.split("_")) if type(i.alleles) is str else []
    return data


def addSequence(cohort_results: CohortMatchResult,
                cohort_answer_seqs: CohortSequences,
                cohort_predit_seqs: CohortSequences) -> None:
    allele_sequence = SeqIO.to_dict(SeqIO.parse("index5/kir_2100_withexon_ab.leftalign.mut01_sequences.fa", "fasta"))
    for allele_name in list(allele_sequence):
        new_allele_name = allele_name
        if new_allele_name[-1].isdigit():
            continue
        if new_allele_name[-1] in ["e", "N"]:
            new_allele_name = new_allele_name[:-1]
        print(allele_name, new_allele_name)
        if new_allele_name not in allele_sequence:
            allele_sequence[new_allele_name] = allele_sequence[allele_name]

    for id, results in cohort_results.items():
        for result in results:
            if result.answer_allele_full:
                result.answer_seq = cohort_answer_seqs[id][result.answer_allele_full].seq
            # if result.answer_allele:
            #     result.answer_seq = allele_sequence[result.answer_allele].seq
            if not result.predit_allele:
                continue

            if id in cohort_predit_seqs and result.predit_allele_full:
                result.predit_seq = cohort_predit_seqs[id][result.predit_allele_full].seq
            else:
                result.predit_seq = allele_sequence[result.predit_allele].seq


def compareCohort(cohort_answer: CohortAlleles,
                  cohort_predit: CohortAlleles,
                  cohort_answer_seqs: CohortSequences = {},
                  cohort_predit_seqs: CohortSequences = {},
                  skip_empty: bool = True,
                  verbose_sample: bool = True,
                  base_compare: bool = False,
                  plot: bool = False) -> dict[str, list[MatchResult]]:
    """
    Compare answer and called alleles from all samples

    Parameters:
      cohort_answer: Reference alleles
      cohort_predit: Predicted alleles
      skip_empty: ignore when there is not such sample id in the predicted cohort
      verbose_sample: Print sample comparison
    """
    samples_id_set = set(cohort_answer.keys())
    if skip_empty:
        samples_id_set = samples_id_set & set(cohort_predit.keys())
    samples_id = sorted(samples_id_set)

    results = {
        sample_id: compareSample(cohort_answer[sample_id],
                                 cohort_predit.get(sample_id, []),
        ) for sample_id in samples_id
    }

    if base_compare:
        addSequence(results, cohort_answer_seqs, cohort_predit_seqs)
        addBaseMatchness(results)

    if verbose_sample:
        print(" === Per sample comparison === ")
        for sample_id, result in results.items():
            print(f"Sample {sample_id}")
            printSummarySample(result)

    if base_compare:
        print(" === Per base comparison === ")
        printSummaryBaseLevel(results)

    print(" === Per gene comparison === ")
    df_gene_acc = printSummaryGeneLevel(results)
    print(" === Per unit comparison === ")
    printSummary(results)
    print(" === Per Resolution comparison === ")
    printSummaryByResolution(results)

    if plot:
        print(" === Plot gene comparison === ")
        showPlot([
            *plotGeneLevelSummary(df_gene_acc),
            *plotGeneLevelSummary(df_gene_acc, order_by_accuracy=True),
        ])
    return results


def compareSample(answer_list: list[str], predict_list: list[str]) -> list[MatchResult]:
    """ Compare two alleles set """
    gene_comparison_result: list[MatchResult] = []

    # Remove output e (In GRAPH-KIR indicate exon-only allele)
    predict_list = [i if not i.endswith("e") else i[:-1] for i in predict_list]

    answer_dict = groupByGene(answer_list)
    predit_dict = groupByGene(predict_list)

    # if KIR2DL5 occurs, treat KIR2DL5A KIR2DL5B in the same gene "KIR2DL5"
    if "KIR2DL5*unresolved" in predict_list:
        answer_dict["KIR2DL5"] = answer_dict.pop("KIR2DL5A", []) + \
                                 answer_dict.pop("KIR2DL5B", [])

    for gene in answer_dict.keys() | predit_dict.keys():
        # if gene in ["KIR2DP1", "KIR3DP1"]:
        #     continue
        gene_comparison_result.extend(compareGene(
            answer_dict[gene],
            predit_dict[gene]
        ))

    return sorted(gene_comparison_result)


def compareGene(a_list: list[str], b_list: list[str]) -> Iterator[MatchResult]:
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
    for allele_b in list(b_list):
        for allele_a in a_list:
            if getAlleleField(allele_a, 7) == getAlleleField(allele_b, 7):
                a_list.remove(allele_a)
                b_list.remove(allele_b)
                yield MatchResult(limitAlleleField(allele_a, 7),  # allele name (remove novel suffix)
                                  limitAlleleField(allele_b, 7),
                                  allele_a, allele_b,
                                  MatchType.MATCH7)
                break

    # Find match 5 digits
    for allele_b in list(b_list):
        for allele_a in a_list:
            if getAlleleField(allele_a, 5) == getAlleleField(allele_b, 5):
                a_list.remove(allele_a)
                b_list.remove(allele_b)
                yield MatchResult(limitAlleleField(allele_a, 7),
                                  limitAlleleField(allele_b, 7),
                                  allele_a, allele_b,
                                  MatchType.MATCH5)
                break

    # Find match 3 digits
    for allele_b in list(b_list):
        for allele_a in a_list:
            if getAlleleField(allele_a, 3) == getAlleleField(allele_b, 3):
                a_list.remove(allele_a)
                b_list.remove(allele_b)
                yield MatchResult(limitAlleleField(allele_a, 7),
                                  limitAlleleField(allele_b, 7),
                                  allele_a, allele_b,
                                  MatchType.MATCH3)
                break

    # Find match < 3 digits
    for allele_a, allele_b in zip(list(a_list), list(b_list)):
        a_list.remove(allele_a)
        b_list.remove(allele_b)
        yield MatchResult(limitAlleleField(allele_a, 7),
                          limitAlleleField(allele_b, 7),
                          allele_a, allele_b,
                          MatchType.MATCHGENE)

    # copy number error
    for allele in a_list:
        yield MatchResult(limitAlleleField(allele, 7), "",
                          allele, "", MatchType.FN)
    for allele in b_list:
        yield MatchResult("", limitAlleleField(allele, 7),
                          "", allele, MatchType.FP)


def printSummarySample(results: list[MatchResult]) -> None:
    """ Print match result """
    for result in results:
        base_compare_str = ""
        if result.answer_allele_length and result.predit_allele_length:
            base_compare_str = f"DIFF {result.base_diff:5.0f} / {result.answer_allele_length:5.0f}"

        if result.match_type == MatchType.MATCH7:
            print(f"{result.answer_allele:18} OK {result.predit_allele:18}", base_compare_str)
        elif result.match_type == MatchType.MATCH5:
            print(f"{result.answer_allele:18} <5 {result.predit_allele:18}", base_compare_str)
        elif result.match_type == MatchType.MATCH3:
            print(f"{result.answer_allele:18} <3 {result.predit_allele:18}", base_compare_str)
        elif result.match_type == MatchType.MATCHGENE:
            print(f"{result.answer_allele:18} <0 {result.predit_allele:18}", base_compare_str)
        elif result.match_type == MatchType.FN:
            print(f"{result.answer_allele:18} <-")
        elif result.match_type == MatchType.FP:
            print(f"{''                  :18} -> {result.predit_allele:18}")


def getBaseMatchness(result: MatchResult) -> dict[str, int]:
    """ subfunction of addBaseMatchness """
    if result.match_type in [
            MatchType.MATCH7,
            MatchType.MATCH5,
            MatchType.MATCH3,
            MatchType.MATCHGENE,
    ]:
        # score = pairwise2.align.globalxx(seq1, seq2, score_only=True)
        score = pairwise2.align.localxx(result.answer_seq, result.predit_seq, score_only=True)

        return {
            'length1': len(result.answer_seq),
            'length2': len(result.predit_seq),
            'score': score,
            'diff': len(result.answer_seq) - score,
        }
    else:
        return {}


def addBaseMatchness(results: CohortMatchResult) -> None:
    """ Concurrently calculate the base difference between alleles """
    with ProcessPoolExecutor() as executor:
        base_results = executor.map(getBaseMatchness, chain.from_iterable(results.values()))
        for result, base_result in zip(chain.from_iterable(results.values()), base_results):
            if not base_result:
                continue
            result.answer_allele_length = base_result["length1"]
            result.predit_allele_length = base_result["length2"]
            result.base_diff            = base_result["diff"]


def printSummary(cohort_result: CohortMatchResult) -> dict[str, int]:
    """
    Summarize the match results in the cohort

    Return:
      summary: The value of each metrics
    """
    summary: dict[str, int] = defaultdict(int)
    for results in cohort_result.values():
        summary['total'] += len([i for i in results if i.match_type != MatchType.FP])
        summary['total_sample'] += 1

        if all(i.match_type == MatchType.FN for i in results):
            summary['fail_allele'] += len(results)
            summary['fail_sample'] += 1
            continue
        for i in results:
            summary[i.match_type.name] += 1

    print(f"Total alleles       = {summary['total']} (sample = {summary['total_sample']})")
    print(f"  * Cannot_called   = {summary['fail_allele']} (sample = {summary['fail_sample']})")
    print(f"  * TP              = {summary[MatchType.MATCH7.name]}")
    print(f"  * Match_5_digits  = {summary[MatchType.MATCH5.name]}")
    print(f"  * Match_3_digits  = {summary[MatchType.MATCH3.name]}")
    print(f"  * Match_gene      = {summary[MatchType.MATCHGENE.name]}")
    print(f"  * Copy_number_err = {summary[MatchType.FN.name] + summary[MatchType.FP.name]}")
    print(f"    * FN = {summary[MatchType.FN.name]}")
    print(f"    * FP = {summary[MatchType.FP.name]}")
    return summary


def calcSummaryByResolution(cohort_results: Iterable[MatchResult]) -> dict[str, int]:
    """ Calculate the correctness by resolution """
    summary = {
        "7digits_total": 0,
        "7digits_correct": 0,
        "5digits_total": 0,
        "5digits_correct": 0,
        "3digits_total": 0,
        "3digits_correct": 0,
        "gene_total": 0,
        "gene_correct": 0,
        "FP": 0,
    }
    for result in cohort_results:
        if not result.answer_allele:
            summary["FP"] += 1
            continue

        if len(getAlleleField(result.answer_allele)) >= 7:
            summary["7digits_total"] += 1
            if result.match_type in [MatchType.MATCH7]:
                summary["7digits_correct"] += 1

        if len(getAlleleField(result.answer_allele)) >= 5:
            summary["5digits_total"] += 1
            if result.match_type in [MatchType.MATCH7, MatchType.MATCH5]:
                summary["5digits_correct"] += 1

        if len(getAlleleField(result.answer_allele)) >= 3:
            summary["3digits_total"] += 1
            if result.match_type in [MatchType.MATCH7, MatchType.MATCH5, MatchType.MATCH3]:
                summary["3digits_correct"] += 1

        if len(getAlleleField(result.answer_allele)) >= 0:
            summary["gene_total"] += 1
            if result.match_type in [MatchType.MATCH7, MatchType.MATCH5, MatchType.MATCH3, MatchType.MATCHGENE]:
                summary["gene_correct"] += 1
    return summary


def printSummaryByResolution(cohort_results: CohortMatchResult) -> dict[str, int]:
    """
    Summarize the match results in the cohort but normalized by answer's resolution

    Return:
      summary: The value of each metrics
    """
    summary = calcSummaryByResolution(result for results in cohort_results.values() for result in results)

    def precisionStr(a: int, b: int) -> str:
        if b != 0:
            return f"{a:5d} / {b:5d} = {a / b:.3f}"
        else:
            return f"{a:5d} / {b:5d} = Nan"

    print(f"7-digits: {precisionStr(summary['7digits_correct'], summary['7digits_total'])}")
    print(f"5-digits: {precisionStr(summary['5digits_correct'], summary['5digits_total'])}")
    print(f"3-digits: {precisionStr(summary['3digits_correct'], summary['3digits_total'])}")
    print(f"Gene:     {precisionStr(summary['gene_correct'],    summary['gene_total'])}")
    print(f"CN:       FP={summary['FP']} FN={summary['gene_total'] - summary['gene_correct']}")
    return summary


def printSummaryBaseLevel(cohort_results: CohortMatchResult) -> pd.DataFrame:
    """
    Summarize the match results in the cohort but normalized by answer's resolution

    Return:
      summary: The value of each metrics
    """
    data = []
    for id, results in cohort_results.items():
        for result in results:
            if not result.answer_allele_length:
                continue
            data.append({
                'diff': result.base_diff,
                'length_ans': result.answer_allele_length,
                'length_pred': result.predit_allele_length,
                'id': id,
            })
            if result.match_type == MatchType.MATCH7:
                data[-1]['match'] = 7
            elif result.match_type == MatchType.MATCH5:
                data[-1]['match'] = 5
            elif result.match_type == MatchType.MATCH3:
                data[-1]['match'] = 3
            else:
                data[-1]['match'] = 0

    df = pd.DataFrame(data)
    print(df.groupby("id").agg({'diff': ["count", "sum", "mean", np.std],
                                'length_ans': ["mean"],
                                'length_pred': ["mean"]}))
    print(df.groupby("match").agg({'diff': ["count", "sum", "mean", np.std],
                                   'length_ans': ["mean"],
                                   'length_pred': ["mean"]}))
    return df


def plotGeneLevelSummary(df: pd.DataFrame,
                         order_by_accuracy: bool = False) -> list[go.Figure]:
    """ Plot gene-level matrics from printSummaryGeneLevel """
    # reorganize the datafrmae
    df_plot = pd.melt(df, id_vars="gene",  # type: ignore
                      value_vars=["7digits_acc", "5digits_acc",
                                  "3digits_acc", "gene_acc"],
                      value_name="accuracy",
                      var_name="level")
    df_plot = df_plot.sort_values(by=["gene", "level"])
    df_plot1 = pd.melt(df, id_vars="gene",  # type: ignore
                       value_vars=["FN", "FP"],
                       value_name="count",
                       var_name="error")
    # order
    color = px.colors.qualitative.Plotly
    df_acc = df_plot.groupby("gene")["accuracy"].mean().reset_index()
    df_acc_sort = df_acc.sort_values(["accuracy", "gene"], ascending=[False, True])
    if order_by_accuracy:
        gene_order = list(df_acc_sort["gene"])
    else:
        gene_order = sorted(df_acc_sort["gene"])
    level_order = ["gene_acc", "3digits_acc", "5digits_acc", "7digits_acc", "FP", "FN"]

    # plot
    fig_accs = px.bar(df_plot, x="gene",  y="accuracy", color="level",
                      color_discrete_sequence=color[:4],
                      category_orders={"level": level_order})
    fig_fnfp = px.bar(df_plot1, x="gene", y="count", color="error", text="count",
                      color_discrete_sequence=color[4:],
                      category_orders={"level": level_order})
    # put them in the same figure
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_traces(list(fig_accs.select_traces()))
    for i in fig_fnfp.select_traces():
        fig.add_trace(i, secondary_y=True)
    fig.update_traces(textposition='outside')
    fig.update_layout(
        xaxis_categoryarray=gene_order,
        xaxis_categoryorder="array",
        yaxis_title="Accuracy",
        yaxis2_title="FN or FP count",
        yaxis2_range=(0, max(df_plot1['count']) * 2),
    )
    return [fig]


def reCalcDigitAcc(df: pd.DataFrame) -> pd.DataFrame:
    df["7digits_acc"] = df["7digits_correct"] / df["7digits_total"]
    df["5digits_acc"] = df["5digits_correct"] / df["5digits_total"]
    df["3digits_acc"] = df["3digits_correct"] / df["3digits_total"]
    df["gene_acc"]    = df["gene_correct"]    / df["gene_total"]
    df["FN"]          = df["gene_total"]      - df["gene_correct"]
    return df


def printSummaryGeneLevel(cohort_results: CohortMatchResult, print_text: bool = True) -> pd.DataFrame:
    """
    printSummaryByResolution but in gene-level mode

    Return:
      summary: The value of each metrics per gene
    """
    # collect result by gene
    result_by_gene = defaultdict(list)
    for results in cohort_results.values():
        for result in results:
            gene = getGeneName(result.answer_allele or result.predit_allele)
            result_by_gene[gene].append(result)

    # per gene accuracy
    summary_by_gene: dict[str, Any] = {}
    for gene in result_by_gene:
        summary_by_gene[gene] = calcSummaryByResolution(result_by_gene[gene])
        summary_by_gene[gene]['gene'] = gene

    # other matrics
    df = pd.DataFrame(summary_by_gene.values())  # type: ignore
    df = reCalcDigitAcc(df)

    if print_text:
        with pd.option_context('display.max_rows', None,
                               'display.max_columns', None,
                               'display.width', 160):  # type: ignore
            print(df.set_index("gene"))
    return df


def compareAlleleWithMethod(cohort_data: dict[str, CohortAlleles]) -> None:
    cohort_results = {}
    for method, dat in cohort_data.items():
        cohort_results[method] = compareCohort(cohort_data["answer"], dat, skip_empty=False, verbose_sample=False)

    cohort_summarys = []
    for method, results in cohort_results.items():
        summary = printSummaryGeneLevel(results, print_text=False)
        # summary = calcSummaryByResolution(result for results in results.values() for result in results)
        summary['method'] = method
        cohort_summarys.append(summary)
    summary = pd.concat(cohort_summarys)

    # print accuracy in differnet perspective
    with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 180):  # type: ignore
        df = summary.groupby(["gene", "method"]).sum()
        print(reCalcDigitAcc(df ).drop(columns=["7digits_total", "5digits_total", "3digits_total", "gene_total"]))

        df1 = summary.groupby(["method"]).sum()
        print(reCalcDigitAcc(df1).drop(columns=["7digits_total", "5digits_total", "3digits_total", "gene_total"]))

        df2 = summary.groupby(["gene"]).sum()
        print(reCalcDigitAcc(df2).drop(columns=["7digits_total", "5digits_total", "3digits_total", "gene_total"]))


if __name__ == "__main__":
    answer = "linnil1_syn/linnil1_syn_s44_summary.csv"
    prefix = "data6/linnil1_syn_s44.{}.30x_s444"
    cohort = [
        {"method": "answer",             "name": f"{answer}"},
        {"method": "ab2dl1s1-pv",        "name": f"{prefix}.index5_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup_assume3DL3.pv.compare_sum.tsv"},
        {"method": "ab2dl1s1-exonfirst", "name": f"{prefix}.index5_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup_assume3DL3.pv_exonfirst_1.2.compare_sum.tsv"},
        {"method": "ping-call",          "name": f"data6/ping_linnil1_syn_{NamePath(prefix).replace_wildcard('_cohort').split('/')[-1]}.result_PING20220527/finalAlleleCalls.merge.tsv"},
        {"method": "gatkir",             "name": f"{NamePath(prefix).replace_wildcard('_merge_called')}.bwa.rg.md.from_linnil1_syn_s44_jg_30x_s444_bwa_rg_md_ploidy_linnil1_syn_s44_same_30x_s444_bwa_rg_md_answer_cn_jg_hc.norm.call_merge.tsv"},
        {"method": "t1k",                "name": f"{NamePath(prefix).replace_wildcard('_mergecall')}.t1k_7_all.tsv"},
        # TODO: exon sequences
    ]

    cohort_data: dict[str, CohortAlleles] = defaultdict(dict)
    for dat in cohort:
        for filename in NamePath(dat['name']).get_input_names():
            print(filename)
            if dat['method'] == "answer":
                cohort_data[dat['method']].update(readAnswerAllele(filename))
            else:
                cohort_data[dat['method']].update(readPredictResult(filename))
    compareAlleleWithMethod(cohort_data)
