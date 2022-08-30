import copy
from pprint import pprint
from typing import TypedDict
from collections import Counter, defaultdict
import numpy as np
import pandas as pd

from kg_eval import readAnswerAllele, getGeneName, CohortAlleles


GeneCN = dict[str, int]


def readCNFile(tsv_file: str) -> GeneCN:
    df = pd.read_csv(tsv_file, sep="\t")
    return dict(zip(map(getGeneName, df['gene']),  df['cn']))


def calCN(alleles: list[str]) -> GeneCN:
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


class CNDiff(TypedDict, total=False):
    gene: str
    total: int
    diff: int
    method: str
    sample_id: str


def compareCN(ans_cn: GeneCN, pred_cn: GeneCN) -> list[CNDiff]:
    """ Compare CN. Return list of {gene, total, diff} in the sample """
    ans_cn = copy.deepcopy(ans_cn)
    pred_cn = copy.deepcopy(pred_cn)

    if "KIR2DL1S1" in pred_cn:
        mergeGene(ans_cn, "KIR2DL1S1", ["KIR2DL1", "KIR2DS1"])
    if "KIR2DL5" in pred_cn:
        mergeGene(ans_cn, "KIR2DL5", ["KIR2DL5A", "KIR2DL5B"])
    if "KIR2DL5AB" in pred_cn:
        mergeGene(ans_cn, "KIR2DL5AB", ["KIR2DL5A", "KIR2DL5B"])
    if "KIR2DS35" in pred_cn:
        mergeGene(ans_cn, "KIR2DS35", ["KIR2DS3", "KIR2DS5"])

    comps = []
    for gene in ans_cn.keys() | pred_cn.keys():
        comp: CNDiff = {'gene': gene, 'total': 0, 'diff': 0}
        if gene in ans_cn:
            comp['total'] += ans_cn[gene]
        comp['diff'] += abs(ans_cn.get(gene, 0) - pred_cn.get(gene, 0))

        if comp['total'] + comp['diff']:  # filter 0
            comps.append(comp)

    return comps


def updateDiff(method: str, sample_id: str, diff: list[CNDiff]) -> list[CNDiff]:
    for i in diff:
        i.update({"method": method, "sample_id": sample_id})
    return diff


if __name__ == "__main__":
    answer = "linnil1_syn_30x_seed87"
    ans = readAnswerAllele(f"{answer}/{answer}.summary.csv")
    data = []
    id_list = list(f"{id:02d}"for id in range(10))
    # for kpi
    # id_list = list(filter(lambda i: i not in [4,5], id_list))
    for id in id_list:
        base = f"data/{answer}.{id}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.linnil1"
        data.extend([{
            'id': id, 'method': "ab2dl1s1_sam_depth",
            'file': f"{base}.cn_sam_depth.tsv",
        },{
            'id': id, 'method': "ab2dl1s1_sam_depth_p75",
            'file': f"{base}.cn_sam_depth_p75.tsv",
        },{
            'id': id, 'method': "ab2dl1s1_sam_depth_p75_kde",
            'file': f"{base}.cn_sam_depth_p75_kde.tsv",
        },{
            'id': id, 'method': "ab2dl1s1_sam_depth_kde_cohert",
            'file': f"{base}.cn_sam_depth_kde_cohert.tsv",
        },{
            'id': id, 'method': "ab2dl1s1_sam_depth_p75_kde_cohert",
            'file': f"{base}.cn_sam_depth_p75_kde_cohert.tsv",
        },{
            'id': id, 'method': "ab2dl1s1_sam_depth_cohert",
            'file': f"{base}.cn_sam_depth_cohert.tsv",
        },{
            'id': id, 'method': "ab2dl1s1_sam_depth_mean",
            'file': f"{base}.cn_sam_depth_cohert.tsv",
        }])

        base = f"data/{answer}.{id}.index_kir_2100_raw.mut01.hisatgenotype.errcorr.linnil1"
        data.extend([{
            'id': id, 'method': "ab_sam_depth",
            'file': f"{base}.cn_sam_depth.tsv",
        },{
            'id': id, 'method': "ab_sam_depth_p75",
            'file': f"{base}.cn_sam_depth_p75.tsv",
        },{
            'id': id, 'method': "ab_sam_depth_cohert",
            'file': f"{base}.cn_sam_depth_cohert.tsv",
        },{
            'id': id, 'method': "ab_sam_depth_kde_cohert",
            'file': f"{base}.cn_sam_depth_kde_cohert.tsv",
        },{
            'id': id, 'method': "ab_sam_depth_p75_kde_cohert",
            'file': f"{base}.cn_sam_depth_p75_kde_cohert.tsv",
        }])

        base = f"data/{answer}.{id}.index_kir_2100_ab.mut01.hisatgenotype.errcorr.linnil1"
        data.extend([{
            'id': id, 'method': "split_sam_depth",
            'file': f"{base}.cn_sam_depth.tsv",
        },{
            'id': id, 'method': "split_sam_depth_p75",
            'file': f"{base}.cn_sam_depth_p75.tsv",
        },{
            'id': id, 'method': "split_sam_depth_cohert",
            'file': f"{base}.cn_sam_depth_cohert.tsv",
        }])

        base = f"data/{answer}_exon.{id}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.linnil1"
        data.extend([{
            'id': id, 'method': "exon-ab2dl1s1_sam_exon_depth",
            'file': f"{base}.cn_sam_exon_depth.tsv",
        },{
            'id': id, 'method': "exon-ab2dl1s1_sam_exon_depth_p75",
            'file': f"{base}.cn_sam_exon_depth_p75.tsv",
        }])
        base = f"data/{answer}.{id}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.linnil1"
        data.extend([{
            'id': id, 'method': "ab2dl1s1_sam_exon_depth",
            'file': f"{base}.cn_sam_exon_depth.tsv",
        },{
            'id': id, 'method': "ab2dl1s1_sam_exon_depth_p75",
            'file': f"{base}.cn_sam_exon_depth_p75.tsv",
        }])


    # print(data)
    results = []
    for i in data:
        comparison = compareCN(calCN(ans[i['id']]), readCNFile(i['file']))
        results.extend(updateDiff(i["method"], i['id'], comparison))

    # PING
    df = pd.read_csv(f"data3/ping_{answer}.result/manualCopyNumberFrame.csv", index_col=0)
    for name, data in df.iterrows():
        # "linnil1_syn_30x_seed87.00.read.",2,1,1,1,0,1,2,2,1,1,1,2,1,1,1,1
        id = name.split('.')[1]
        comparison = compareCN(calCN(ans[id]), dict(data))
        results.extend(updateDiff("PING2", id, comparison))

    # KPI
    df = pd.read_csv(f"data3/{answer}_merge_cn.kpi_prediction.csv", index_col=0)
    for name, data in df.iterrows():
        id = name.split('.')[1]
        comparison = compareCN(calCN(ans[id]), dict(data))
        results.extend(updateDiff("KPI", id, comparison))

    # GATKIR
    df = pd.read_csv(f"data3/{answer}_merge_depth.bwa.rg.md.coverage.depth_per_gene.ploidy.csv", index_col=0).T
    for name, data in df.iterrows():
        id = name
        comparison = compareCN(calCN(ans[id]), dict(data))
        results.extend(updateDiff("GATKIR", id, comparison))

    # summary
    df = pd.DataFrame(results)
    df = df.loc[df['sample_id'].isin(id_list)]
    # df = df.loc[np.logical_not(df['gene'].str.contains("KIR3DL2|KIR3DP1"))]
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # type: ignore
        df1 = df.groupby(["method", "gene"]).sum()
        df1['acc'] = 1 - df1['diff'] / df1['total']
        print(df1)

    df2 = df.groupby(["method"]).sum()
    df2['acc'] = 1 - df2['diff'] / df2['total']
    print(df2)
