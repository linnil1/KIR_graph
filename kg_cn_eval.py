import pandas as pd
from kg_eval import EvaluateKIR, getGeneName
from collections import Counter, defaultdict
import copy
from pprint import pprint


def readCNFile(f):
    df = pd.read_csv(f, sep="\t")
    return dict(zip(map(getGeneName, df['gene']),  df['cn']))


def mergeGene(gene_cn, gene_to, gene_froms):
    cn = 0
    for gene in gene_froms:
        if gene in gene_cn:
            cn += gene_cn[gene]
            gene_cn.pop(gene)
    gene_cn[gene_to] = cn
    return gene_cn


def compareCN(ans_cn, pred_cn):
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
        comp = {'gene': gene, 'total': 0, 'diff': 0}
        if gene in ans_cn:
            comp['total'] += ans_cn[gene]
        comp['diff'] += abs(ans_cn.get(gene, 0) - pred_cn.get(gene, 0))

        if comp['total'] + comp['diff']:  # filter 0
            comps.append(comp)

    return comps


def mergeDict(d_from, d_to):
    for i, cn in d_from.items():
        d_to[i] += cn


answer = "linnil1_syn_30x_seed87"
ans = EvaluateKIR(f"{answer}/{answer}.summary.csv")
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



print(data)
method = defaultdict(list)
for i in data:
    comparison = compareCN(ans.getAnsCN(i['id']), readCNFile(i['file']))
    method[i['method']].append(comparison)

# PING
df = pd.read_csv(f"data3/ping_{answer}.result/manualCopyNumberFrame.csv", index_col=0)
for name, data in df.iterrows():
    # "linnil1_syn_30x_seed87.00.read.",2,1,1,1,0,1,2,2,1,1,1,2,1,1,1,1
    id = name.split('.')[1]
    if id not in id_list:
        continue
    comparison = compareCN(ans.getAnsCN(id), dict(data))
    method["PING2"].append(comparison)

# KPI
df = pd.read_csv(f"data3/{answer}_merge_cn.kpi_prediction.csv", index_col=0)
for name, data in df.iterrows():
    id = name.split('.')[1]
    if id not in id_list:
        continue
    comparison = compareCN(ans.getAnsCN(id), dict(data))
    method["KPI"].append(comparison)

# GATKIR
df = pd.read_csv(f"data3/{answer}_merge_depth.bwa.rg.md.coverage.depth_per_gene.ploidy.csv", index_col=0).T
for name, data in df.iterrows():
    id = name
    if id not in id_list:
        continue
    comparison = compareCN(ans.getAnsCN(id), dict(data))
    method["GATKIR"].append(comparison)


merged_method = []
for i, compare_list in method.items():
    df = pd.DataFrame([j for i in compare_list for j in i])
    df = df.groupby("gene").sum()
    print(i)
    print(df)
    # df = df.drop(index=["KIR3DL2", "KIR3DP1"])
    df['acc'] = 1 - df['diff'] / df['total']
    summary = {
        'method': i,
        'diff': df['diff'].sum(),
        'total': df['total'].sum(),
        'acc': 1 - df['diff'].sum() / df['total'].sum(),
    }
    merged_method.append(summary)

print(pd.DataFrame(merged_method))
