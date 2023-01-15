import copy
from pprint import pprint
from typing import TypedDict
from collections import Counter, defaultdict
import numpy as np
import pandas as pd

from namepipe import NamePath
from graphkir.utils import getGeneName
from kg_eval import readAnswerAllele, CohortAlleles, readPredictResult
from kg_utils import getAnswerFile


GeneCN = dict[str, int]
CohortGeneCN = dict[str, GeneCN]


class CNDiff(TypedDict, total=False):
    """The dict before insert into Dataframe"""

    gene: str
    total: int
    diff: int
    diff_abs: int
    method: str
    sample_id: str


def readCNFile(tsv_file: str) -> GeneCN:
    """Read CN file"""
    df = pd.read_csv(tsv_file, sep="\t")
    return dict(zip(map(getGeneName, df["gene"]), df["cn"]))


def allele2CN(alleles: list[str]) -> GeneCN:
    """list of alleles -> CN of each gene"""
    return Counter(map(getGeneName, alleles))


def mergeGene(gene_cn: GeneCN, gene_to: str, gene_froms: list[str]) -> GeneCN:
    """
    Treat all gene in `gene_froms` as same gene `gene_to`.

    Example:
      Input: KIR2DL1 cn=1, KIR2DS1 cn=2
      Output: KIR2DL1S1 cn=3
    """
    cn = 0
    for gene in gene_froms:
        if gene in gene_cn:
            cn += gene_cn[gene]
            gene_cn.pop(gene)
    gene_cn[gene_to] = cn
    return gene_cn


def compareCN(ans_cn: GeneCN, pred_cn: GeneCN) -> list[CNDiff]:
    """Compare CN. Return list of {gene, total, diff} in the sample"""
    ans_cn = copy.deepcopy(ans_cn)
    pred_cn = copy.deepcopy(pred_cn)

    if "KIR2DL1S1" in pred_cn:
        mergeGene(ans_cn, "KIR2DL1S1",         ["KIR2DL1", "KIR2DS1"])
    if "KIR2DL5" in pred_cn:
        mergeGene(ans_cn, "KIR2DL5",           ["KIR2DL5A", "KIR2DL5B"])
    if "KIR2DL5AB" in pred_cn:
        mergeGene(ans_cn, "KIR2DL5AB",         ["KIR2DL5A", "KIR2DL5B"])
    if "KIR2DS35" in pred_cn:
        mergeGene(ans_cn, "KIR2DS35",          ["KIR2DS3", "KIR2DS5"])
    if "KIR2DL5A;KIR2DL5B" in pred_cn:
        mergeGene(ans_cn, "KIR2DL5A;KIR2DL5B", ["KIR2DL5A", "KIR2DL5B"])
    if "KIR2DS3;KIR2DS5" in pred_cn:
        mergeGene(ans_cn, "KIR2DS3;KIR2DS5",   ["KIR2DS3", "KIR2DS5"])

    comps = []
    for gene in ans_cn.keys() | pred_cn.keys():
        comp: CNDiff = {"gene": gene, "total": 0, "diff": 0, "diff_abs": 0}
        if gene in ans_cn:
            comp["total"] += ans_cn[gene]
        diff = int(ans_cn.get(gene, 0)) - int(pred_cn.get(gene, 0))
        comp["diff"] += diff
        comp["diff_abs"] += abs(diff)

        if comp["total"] + comp["diff_abs"]:  # filter 0
            comps.append(comp)

    return comps


def compareCNCohort(
    cohort_cn_ans: CohortGeneCN, cohort_cn_pred: CohortGeneCN, method: str
) -> list[CNDiff]:
    """The difference of answer and prediction"""
    diffs_cohort = []
    for sample_id, ans_cn in cohort_cn_ans.items():
        diffs = compareCN(ans_cn, cohort_cn_pred.get(sample_id, {}))
        for diff in diffs:
            diff["sample_id"] = sample_id
            diff["method"] = method
        diffs_cohort.extend(diffs)
    return diffs_cohort


def compareCNwithMethods(
    method_cohort_cn: dict[str, CohortGeneCN],
    print_per_method_gene: bool = True,
    print_per_method_sample: bool = True,
    print_per_method: bool = True,
) -> None:
    """Compare the CN per sample, per gene, per method"""
    # from pprint import pprint
    # pprint(cohort_data)
    results = []
    for method, result in method_cohort_cn.items():
        results.extend(compareCNCohort(method_cohort_cn["answer"], result, method))
    # pprint(results)

    # summary
    df = pd.DataFrame(results)
    # remove samples
    # df = df.loc[np.logical_not(df['sample_id'].isin(["09"]))]
    # remove samples (for HPRC data)
    # df = df.loc[np.logical_not(df['sample_id'].isin(["HG00733", "HG01123", "HG02486", "HG02559", "NA19240"]))]
    # remove genes  (GATKIR)
    # df = df.loc[np.logical_not(df['gene'].str.contains("KIR3DL2|KIR3DP1"))]
    with pd.option_context("display.max_rows", None, "display.max_columns", None):  # type: ignore
        if print_per_method_gene:
            df1 = df.groupby(["method", "gene"]).sum()
            df1["acc"] = 1 - df1["diff_abs"] / df1["total"]
            print(df1)

        if print_per_method_sample:
            df3 = df.groupby(["method", "sample_id"]).sum()
            print(df3)

    if print_per_method:
        df2 = df.groupby(["method"]).sum()
        df2["acc"] = 1 - df2["diff_abs"] / df2["total"]
        print(df2)


def readAnyMethodToCN(filename: NamePath | str, method: str) -> CohortGeneCN:
    """CN reading factory"""
    print(method, filename)
    if method == "answer":
        return {
            i: allele2CN(j)
            for i, j in readAnswerAllele(getAnswerFile(filename)).items()
        }

    if method == "allele":  # linnil1's allele format
        cohort = {}
        for f in NamePath(filename).get_input_names():
            cohort.update({
                i: allele2CN(j) for i, j in readPredictResult(filename).items()
            })
        return cohort

    if method == "cn":  # linnil1's CN
        cohort = {}
        for f in NamePath(filename).get_input_names():
            if not f.endswith(".tsv"):
                f += ".tsv"
            cohort[f.template_args[0]] = readCNFile(f)
        return cohort

    if method == "ping":  # ping's intermediate CN result
        df = pd.read_csv(filename, index_col=0)
        cohort = {}
        for name, data in df.iterrows():
            # "linnil1_syn_30x_seed87.00.read.",2,1,1,1,0,1,2,2,1,1,1,2,1,1,1,1
            id = str(name).split(".")[1]  # TOOD: temporary works
            cohort[id] = dict(data)
        return cohort

    if method == "sakauekir":
        df = pd.read_csv(filename, index_col=0).T
        cohort = {}
        for name, data in df.iterrows():
            cohort[str(name)] = dict(data)
        return cohort

    if method == "kpi":
        df = pd.read_csv(filename, index_col=0).T
        cohort = {}
        for name, data in df.iterrows():
            id = str(name)  # .split('.')[1]
            cohort[id] = dict(data)
        return cohort

    raise ValueError(f"Method {method} not found")


def compareCNResult(
    input_name: NamePath, sample_name: NamePath, method: str = "cn"
) -> NamePath:
    """Print comparison in CN"""
    compareCNwithMethods(
        {
            "answer": readAnyMethodToCN(sample_name, "answer"),
            "linnil1": readAnyMethodToCN(input_name, "cn"),
        }
    )
    return input_name


if __name__ == "__main__":
    answer = "linnil1_syn/linnil1_syn_s44_summary"
    prefix = "data/linnil1_syn_s44.{}.30x_s444"
    answer = "linnil1_syn/linnil1_syn_s2022_summary"
    prefix = "data/linnil1_syn_s2022.{}.30x_s1031"
    real = False
    if not real:
        cohort = [
            {"method": "answer",
             "name": f"{answer}"},
            # {"method": "ab2dl1s1-dev008-call",
            #  "name": f"{prefix}.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup_assume3DL3.pv.compare_sum.tsv"},
            # {"method": "ab2dl1s1",
            #  "name": f"{prefix}.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup_dev0.06_assume3DL3.tsv"},
            # {"method": "ab2dl1s1-sg",
            #  "name": f"{prefix}.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup_smallg_dev0.06_assume3DL3.tsv"},
            {"method": "ab2dl1s1-b2",
             "name": f"{prefix}.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup_b2_assume3DL3.tsv"},
            {"method": "ab2dl1s1-b2-call",
             "name": f"{prefix}.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup_b2_assume3DL3.pv.compare_sum.var_errcorr.top600.tsv"},
            {"method": "ab2dl1s1-dev6",
             "name": f"{prefix}.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup_dev0.06_assume3DL3.tsv"},
            {"method": "ab2dl1s1-dev6-call",
             "name": f"{prefix}.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup_dev0.06_assume3DL3.pv.compare_sum.var_errcorr.top600.tsv"},
            # {"method": "ab2dl1s1-exon-call", "name": f"{prefix}.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup_dev0.06_assume3DL3.pv_exonfirst_1.compare_sum.var_errcorr.top600.tsv"},
            # {"method": "ab2dl1s1-cohort", "name": f"{prefix}.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup_dev0.06.cohort.tsv"},
            # {"method": "ab2dl1s1-cohort-dev008", "name": f"{prefix}.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup.cohort.tsv"},
            # {"method": "ab2dl1s1-cohortsg", "name": f"{prefix}.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup_smallg_dev0.06.cohort.tsv"},
            # {"method": "ab2dl1s1-median", "name": f"{prefix}.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.median.CNgroup_assume3DL3.tsv"},
            # {"method": "ab2dl1s1-kde",    "name": f"{prefix}.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.kde.cohort.tsv"},
            # {"method": "ab-cohort",       "name": f"{prefix}.index_kir_2100_withexon_ab.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup.cohort.tsv"},
            # {"method": "ab",              "name": f"{prefix}.index_kir_2100_withexon_ab.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup_assume3DL3.tsv"},
            # {"method": "split-cohort",    "name": f"{prefix}.index_kir_2100_withexon_split.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup.cohort.tsv"},
            # {"method": "split",           "name": f"{prefix}.index_kir_2100_withexon_split.leftalign.mut01.graph.variant.noerrcorr.no_multi.depth.p75.CNgroup_assume3DL3.tsv"},
            {"method": "ping",          "name": f"{NamePath(prefix).replace_wildcard('_pingsample')}.result_ping_20220527/manualCopyNumberFrame.csv", },
            {"method": "ping-call",     "name": f"{NamePath(prefix).replace_wildcard('_pingsample')}.result_ping_20220527.merge.tsv", },
            {"method": "ping-wgs",      "name": f"{NamePath(prefix).replace_wildcard('_pingsample')}.result_ping_wgs/manualCopyNumberFrame.csv", },
            {"method": "ping-wgs-call", "name": f"{NamePath(prefix).replace_wildcard('_pingsample')}.result_ping_wgs.merge.tsv", },
            {"method": "sakauekir",     "name": f"{NamePath(prefix).replace_wildcard('_merge_depth')}.sakauekir_v1_0_0.bwa.rg.md.coverage.depth.ploidy.csv", },
            {"method": "kpi",           "name": f"{NamePath(prefix).replace_wildcard('_merge_cn')}.kpi_kpi_v1_1_1_prediction.csv", },
            {"method": "t1k-call",      "name": f"{NamePath(prefix).replace_wildcard('_mergecall')}.t1k_t1k_v1_0_1_ipd_2100.dig7.tsv", },
            # TODO: exon-only
        ]
    else:
        answer = "hprc_summary"
        prefix = "data_real/hprc.{}.index_hs37d5.bwa.part_strict"  # .annot_read"
        prefix_graph = ".index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim"
        cohort = [
            {"method": "answer", "name": f"{answer}"},
            {"method": "ab2dl1s1-dev008",
             "name": f"{prefix}{prefix_graph}.variant.noerrcorr.no_multi.depth.p75.CNgroup_assume3DL3.tsv"},
            {"method": "ab2dl1s1-b2",
             "name": f"{prefix}{prefix_graph}.variant.noerrcorr.no_multi.depth.p75.CNgroup_b2_assume3DL3.tsv"},
            {"method": "ab2dl1s1-b2-call",
             "name": f"{prefix}{prefix_graph}.variant.noerrcorr.no_multi.depth.p75.CNgroup_b2_assume3DL3.pv.compare_sum.var_errcorr.top600.tsv"},
            {"method": "ping",          "name": f"{NamePath(prefix).replace_wildcard('_pingsample')}.result_ping_20220527/manualCopyNumberFrame.csv"},
            {"method": "ping-call",     "name": f"{NamePath(prefix).replace_wildcard('_pingsample')}.result_ping_20220527.merge.tsv"},
            {"method": "ping-wgs",      "name": f"{NamePath(prefix).replace_wildcard('_pingsample')}.result_ping_wgs/manualCopyNumberFrame.csv"},
            {"method": "ping-wgs-call", "name": f"{NamePath(prefix).replace_wildcard('_pingsample')}.result_ping_wgs.merge.tsv"},
        ]

    cohort_data: dict[str, CohortGeneCN] = defaultdict(dict)
    for dat in cohort:
        print(dat)
        short_method_name = dat["method"]
        if short_method_name.endswith("-call"):
            short_method_name = "allele"
        if "ping" in short_method_name:
            short_method_name = "ping"
        elif short_method_name not in ["sakauekir", "kpi", "ping", "answer", "allele"]:
            short_method_name = "cn"
        cohort_data[dat["method"]] = readAnyMethodToCN(dat["name"], short_method_name)

    compareCNwithMethods(cohort_data)
