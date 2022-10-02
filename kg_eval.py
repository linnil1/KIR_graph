import re
from enum import Enum
from typing import Iterator, Callable, Iterable, Any
from dataclasses import dataclass
from collections import defaultdict
import pandas as pd

from graphkir.utils import getGeneName, getAlleleField


CohortAlleles = dict[str, list[str]]  # it's ordered dict too


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
    answer_allele: str = ""  # order is important
    predit_allele: str = ""
    match_type: MatchType = MatchType.NONE

    def __lt__(self, other: 'MatchResult') -> bool:
        return (self.answer_allele or self.predit_allele) \
               < (other.answer_allele or other.predit_allele)


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
    data = pd.read_csv(tsv_file, sep='\t')
    return {extract_func(i.name): sorted(i.alleles.split("_"))
                for i in data.itertuples()}


def compareCohort(cohort_answer: CohortAlleles, cohort_predit: CohortAlleles,
                  skip_empty: bool = True, verbose_sample: bool = True) -> None:
    """
    Compare answer and called alleles from all samples

    Parameters:
      cohort_answer: Reference alleles
      cohort_predit: Predicted alleles
      skip_empty: ignore when there is not such sample id in the predicted cohort
      verbose_sample: Print sample comparison
    """
    samples_id = set(cohort_answer.keys())
    if skip_empty:
        samples_id = samples_id & set(cohort_predit.keys())

    results = []
    for sample_id in sorted(samples_id):
        results.append(compareSample(
            cohort_answer[sample_id],
            cohort_predit.get(sample_id, []),
        ))
        if verbose_sample:
            print(f"Sample {sample_id}")
            printSummarySample(results[-1])

    printSummary(results)
    printSummaryByResolution(results)
    printSummaryGeneLevel(results)


def compareSample(answer_list: list[str], predict_list: list[str]) -> list[MatchResult]:
    """ Compare two alleles set """
    gene_comparison_result: list[MatchResult] = []
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
    for allele in list(b_list):
        if allele in a_list:
            a_list.remove(allele)
            b_list.remove(allele)
            yield MatchResult(allele, allele, MatchType.MATCH7)

    # Find match 5 digits
    for allele_b in list(b_list):
        for allele_a in a_list:
            if getAlleleField(allele_a, 5) == getAlleleField(allele_b, 5):
                a_list.remove(allele_a)
                b_list.remove(allele_b)
                yield MatchResult(allele_a, allele_b, MatchType.MATCH5)
                break

    # Find match 3 digits
    for allele_b in list(b_list):
        for allele_a in a_list:
            if getAlleleField(allele_a, 3) == getAlleleField(allele_b, 3):
                a_list.remove(allele_a)
                b_list.remove(allele_b)
                yield MatchResult(allele_a, allele_b, MatchType.MATCH3)
                break

    # Find match < 3 digits
    for allele_a, allele_b in zip(list(a_list), list(b_list)):
        a_list.remove(allele_a)
        b_list.remove(allele_b)
        yield MatchResult(allele_a, allele_b, MatchType.MATCHGENE)

    # copy number error
    for allele in a_list:
        yield MatchResult(allele, "", MatchType.FN)
    for allele in b_list:
        yield MatchResult("", allele, MatchType.FP)


def printSummarySample(results: list[MatchResult]) -> None:
    """ Print match result """
    for result in results:
        if result.match_type == MatchType.MATCH7:
            print(f"{result.answer_allele:18} OK {result.predit_allele:18}")
        elif result.match_type == MatchType.MATCH5:
            print(f"{result.answer_allele:18} <5 {result.predit_allele:18}")
        elif result.match_type == MatchType.MATCH3:
            print(f"{result.answer_allele:18} <3 {result.predit_allele:18}")
        elif result.match_type == MatchType.MATCHGENE:
            print(f"{result.answer_allele:18} <0 {result.predit_allele:18}")
        elif result.match_type == MatchType.FN:
            print(f"{result.answer_allele:18} <-")
        elif result.match_type == MatchType.FP:
            print(f"{''                  :18} -> {result.predit_allele:18}")


def printSummary(cohort_result: list[list[MatchResult]]) -> dict[str, int]:
    """
    Summarize the match results in the cohort

    Return:
      summary: The value of each metrics
    """
    summary: dict[str, int] = defaultdict(int)
    for results in cohort_result:
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

        if len(getAlleleField(result.answer_allele)) == 7:
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



def printSummaryByResolution(cohort_results: list[list[MatchResult]]) -> dict[str, int]:
    """
    Summarize the match results in the cohort but normalized by answer's resolution

    Return:
      summary: The value of each metrics
    """
    summary = calcSummaryByResolution(result for results in cohort_results for result in results)

    def precisionStr(a: int, b: int) -> str:
        return f"{a:5d} / {b:5d} = {a / b:.3f}"

    print(f"7-digits: {precisionStr(summary['7digits_correct'], summary['7digits_total'])}")
    print(f"5-digits: {precisionStr(summary['5digits_correct'], summary['5digits_total'])}")
    print(f"3-digits: {precisionStr(summary['3digits_correct'], summary['3digits_total'])}")
    print(f"Gene:     {precisionStr(summary['gene_correct'],    summary['gene_total'])}")
    print(f"CN:       FP={summary['FP']} FN={summary['gene_total'] - summary['gene_correct']}")
    return summary


def plotGeneLevelSummary(df: pd.DataFrame) -> None:
    """ Plot gene-level matrics from printSummaryGeneLevel """
    from graphkir.plot import showPlot
    from plotly.subplots import make_subplots
    import plotly.express as px
    # reorganize the datafrmae
    df_plot = pd.melt(df, id_vars="gene",  # type: ignore
                      value_vars=['7digits_acc', "5digits_acc",
                                  "3digits_acc", "gene_acc"],
                      value_name="accuracy",
                      var_name="level")
    df_plot = df_plot.sort_values(by=["gene", "level"])
    df_plot1 = pd.melt(df, id_vars="gene",  # type: ignore
                       value_vars=["FN", "FP"],
                       value_name="count",
                       var_name="error")
    # plot in the same figure
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_traces(list(px.line(df_plot, x="gene", y="accuracy", color="level").select_traces()))
    for i in px.bar(df_plot1, x="gene", y="count", color="error", opacity=0.5).select_traces():
        fig.add_trace(i, secondary_y=True)
    fig.update_layout(
        yaxis_title="Accuracy",
        yaxis2_title="FN or FP count",
        yaxis2_range=(0, max(df_plot1['count']) * 2),
    )
    showPlot([fig])


def printSummaryGeneLevel(cohort_results: list[list[MatchResult]]) -> dict[str, dict[str, int]]:
    """
    printSummaryByResolution but in gene-level mode

    Return:
      summary: The value of each metrics per gene
    """
    # collect result by gene
    result_by_gene = defaultdict(list)
    for results in cohort_results:
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
    df["7digits_acc"] = df["7digits_correct"] / df["7digits_total"]
    df["5digits_acc"] = df["5digits_correct"] / df["5digits_total"]
    df["3digits_acc"] = df["3digits_correct"] / df["3digits_total"]
    df["gene_acc"]    = df["gene_correct"]    / df["gene_total"]
    df["FN"]          = df["gene_total"]      - df["gene_correct"]

    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # type: ignore
    #     print(df)
    # plotGeneLevelSummary(df)
    return summary_by_gene


if __name__ == "__main__":
    answer = "linnil1_syn_30x_seed87"

    data = [{
        'type': 'ping',
        'file': f"data3/ping_{answer}.result/finalAlleleCalls.csv",
    }, {
        'type': 'ping-20220527',
        'file': f"data3/ping_{answer}.resultPING20220527/finalAlleleCalls.csv",
    }, {
        'type': 'hisat-271-ab',
        'file': f"data/{answer}_merge.index_kir_271_raw.mut01.hisatgenotype.errcorr.linnil1.cn_sam_depth_p75.type_likelihood.tsv"
    }, {
       'type': 'hisat-271-ab-2dl1s1',
       'file': f"data/{answer}_merge.index_kir_271_2dl1s1.mut01.hisatgenotype.errcorr.linnil1.cn_sam_depth_p75.type_likelihood.tsv"
    }, {
        'type': 'hisat-290-ab',
        'file': f"data/{answer}_merge.index_kir_290_raw.mut01.hisatgenotype.errcorr.linnil1.cn_sam_depth_p75.type_likelihood.tsv"
    }, {
       'type': 'hisat-290-ab-2dl1s1',
       'file': f"data/{answer}_merge.index_kir_290_2dl1s1.mut01.hisatgenotype.errcorr.linnil1.cn_sam_depth_p75.type_likelihood.tsv"
    }, {
       'type': 'GATKIR',
       'file': f"data3/{answer}_merge_called.bwa.rg.md.from_linnil1_syn_30x_seed87_jg_bwa_rg_md_ploidy_linnil1_syn_30x_seed87_answer_cn_jg_hc.norm.call_merge.tsv",
    }, {
       'type': 'GATKIR-all',
       'file': f"data3/{answer}_merge_called_full.bwa.rg.md.from_linnil1_syn_30x_seed87_jg_bwa_rg_md_ploidy_linnil1_syn_30x_seed87_answer_cn_jg_hc.norm.call_merge.tsv",
    }]

    data.extend([{
        'type': 'hisat-ab',
        'file': f"data/{answer}_merge.index_kir_2100_raw.mut01.hisatgenotype.errcorr.linnil1.cn_sam_depth_p75.type_likelihood.tsv"
    }, {
        'type': 'hisat-ab-2dl1s1',
        'file': f"data/{answer}_merge.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.linnil1.cn_sam_depth_p75.type_likelihood.tsv"
    }, {
        'type': 'hisat',
        'file': f"data/{answer}_merge.index_kir_2100_ab.mut01.hisatgenotype.errcorr.linnil1.cn_sam_depth_p75.type_likelihood.tsv"
    }, {
        'type': 'hisat-ab-2dl1s1-report',
        'file': f"data/{answer}_merge.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.linnil1.cn_sam_depth_p75.type_hisat.tsv"
    }, {
        'type': 'hisat-ab-2dl1s1-multi',
        'file': f"data/{answer}_merge.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.linnil1.cn_sam_depth_p75.type_likelihood_multi.tsv"
    }, {
        'type': 'hisat-ab-2dl1s1-no-errorcorr',
        'file': f"data/{answer}_merge.index_kir_2100_2dl1s1.mut01.hisatgenotype.linnil1.cn_sam_depth_p75.type_likelihood.tsv"
    }, {
        'type': 'hisat-ab-2dl1s1',
        'file': f"data/{answer}_merge.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.linnil1.cn_sam_depth_p75.type_likelihood.tsv"
    }])

    answer += "_exon"
    data.extend([{
        'type': 'ping-exon',
        'file': f"data3/ping_{answer}.result/finalAlleleCalls.csv",
    }, {
        'type': 'ping-20220527-exon',
        'file': f"data3/ping_{answer}.resultPING20220527/finalAlleleCalls.csv",
    }, {
        'type': 'hisat-ab-exon',
        'file': f"data/{answer}_merge.index_kir_2100_raw.mut01.hisatgenotype.errcorr.linnil1.cn_sam_exon_depth_p75.type_likelihood.tsv"
    }, {
        'type': 'hisat-ab-2dl1s1-exon',
        'file': f"data/{answer}_merge.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.linnil1.cn_sam_exon_depth_p75.type_likelihood.tsv"
    # }, {
    #     'type': 'hisat-exon',
    #     'file': f"data/{answer}_merge.index_kir_2100_ab.mut01.hisatgenotype.errcorr.linnil1.cn_sam_exon_depth_p75.type_likelihood.tsv"
    }])

    from ping import readPingResult
    cohort = readAnswerAllele(f"{answer}/{answer}.summary.csv")

    for dat in data:
        print(dat['type'])
        if "ping" in dat['type']:
            called = readPingResult(dat['file'])
        else:
            called = readPredictResult(dat['file'])
        compareCohort(cohort, called, skip_empty=False, verbose_sample=False)
