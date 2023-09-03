# #: identifying an allele having nonsynonymous variant(s) in CDS, segmental deletion or fusion with another gene
# =: identifying an allele having synonymous variant(s) in CDS
# $: identifying an allele having variant(s) in non-CDS
# +: identifying a genomic allele matching an IPD-KIR CDS-only allele
# ~: a separator of two genes, usually upstream/downstream or parts of duplication genes
# -: a connector of fusion genes
# _e: identifying only designated exons of an allele
from typing import Any
from pprint import pprint
from functools import partial
import pandas as pd
from kg_eval import readPredictResult


def fusionRewrite(alleles):
    # KIR2DS4*00101_e124567-3DL1*03501_e8
    for allele in alleles:
        if "_e" not in allele:
            yield allele
            continue
        # print(allele)
        allele_rewrite = "e".join(
            [allele.split("_e")[0] for allele in allele.split("-")]
        )
        # print(allele_rewrite)
        yield allele_rewrite


def novelExonRewrite(alleles):
    # KIR3DP1*004(99.89%)
    # KIR3DL3# (99.92% of 00902)
    for allele in alleles:
        if "(" not in allele:
            yield allele
            continue
        # print(allele)
        if "#" in allele:
            allele_rewrite = (
                allele.split("#")[0]
                + "*"
                + allele.split("of")[1].split(")")[0].strip()
                + "#"
            )
        else:
            allele_rewrite = allele.split("(")[0] + "#"
        # print(allele_rewrite)
        yield allele_rewrite


def mergeAlleleField(sample: pd.Series) -> pd.Series:
    allele_fields = sample.iloc[3:]
    allele_fields = allele_fields[allele_fields.notna()]
    # print("Before", allele_fields)
    allele_fields = allele_fields.str.findall(r"KIR[^~]*").explode()
    allele_fields = allele_fields[allele_fields.notna()]
    allele_fields = fusionRewrite(allele_fields)
    allele_fields = novelExonRewrite(allele_fields)
    allele_fields = list(allele_fields)
    # print("After\n", "\n".join(allele_fields), sep="")
    return pd.Series([allele_fields], dtype="object")


def mergeHaplo(sample: pd.DataFrame) -> pd.Series:
    # sort by KIR name, but the order will not preserved
    sample["alleles"] = sample["alleles"].apply(sorted)
    return pd.Series(
        {
            "haplos": "+".join(sample["alleles"].str.len().astype(str)),
            "alleles": "_".join(sample["alleles"].sum()),
        }
    )


def saveKelvinHPRC():
    hprc_kir_df = pd.read_csv(input_csv, header=None)
    hprc_kir_df = hprc_kir_df.drop(index=[0, 2, 3, 4, 5])  # header & hg002 extra
    hprc_kir_df.loc[1, 0] = "hg002"
    hprc_kir_df["alleles"] = hprc_kir_df.apply(mergeAlleleField, axis=1)
    hprc_kir_df = hprc_kir_df[hprc_kir_df["alleles"].str.len() > 0]
    hprc_kir_df = hprc_kir_df.rename(columns={0: "id", 1: "PM"})
    hprc_kir_df = hprc_kir_df[["id", "PM", "alleles"]]
    print(hprc_kir_df)

    # transform to our format
    hprc_kir_df_sample = hprc_kir_df.groupby("id", as_index=False).apply(mergeHaplo)
    hprc_kir_df_sample["name"] = hprc_kir_df_sample["id"]
    hprc_kir_df_sample["id"] = hprc_kir_df_sample["id"].str.upper()
    hprc_kir_df_sample["name"] = hprc_kir_df_sample["name"].str.upper()
    print(hprc_kir_df_sample)
    hprc_kir_df_sample.to_csv(output_csv, sep="\t", index=False)


def removeNovel(df, novel_res=-1):
    genes = set()
    if novel_res >= 0:
        for allele in set(df[df["allele1"].str.contains("e", na=False, regex=False)]["allele1"]):
            genes.update(i.split("*")[0] for i in allele.split("e"))
    if novel_res >= 3:
        genes.update(df[df["allele1"].str.contains("#", na=False, regex=False)]["gene"])
    if novel_res >= 5:
        genes.update(df[df["allele1"].str.contains("=", na=False, regex=False)]["gene"])
        genes.update(df[df["allele1"].str.contains("+", na=False, regex=False)]["gene"])
    if novel_res >= 7:
        genes.update(df[df["allele1"].str.contains("$", na=False, regex=False)]["gene"])
    df = df[~df["gene"].isin(genes)]
    return df

from star_allele_comp import compare_method, table_summarize, compact_summary

if __name__ == "__main__":
    input_csv = "kelvin_kir_v1_0.csv"
    output_csv = "hprc_summary_v1_0_e.tsv"
    # saveKelvinHPRC()
    methods = [
        {
            "name": "answer",
            "file": output_csv,
        },
        {
            "name": "graphkir_full",
            "file": "data_real/hprc_merge.index_hs37d5.bwa.part_strict.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim.variant.noerrcorr.no_multi.depth.p75.CNgroup_b2_assume3DL3.pv.compare_sum.var_errcorr.top600.tsv",
        },
        {
            "name": "graphkir_exonfirst",
            "file": "data_real/hprc_merge.index_hs37d5.bwa.part_strict.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim.variant.noerrcorr.no_multi.depth.p75.CNgroup_b2_assume3DL3.pv_exonfirst_1.compare_sum.var_errcorr.top600.tsv",
        },
        {
            "name": "ping",
            "file": "data_real/hprc_pingsample.index_hs37d5.bwa.part_strict.result_ping_20220527.merge.tsv",
        },
        {
            "name": "ping-wgs",
            "file": "data_real/hprc_pingsample.index_hs37d5.bwa.part_strict.result_ping_wgs.merge.tsv",
        },
    ]

    cohort = {}
    for method in methods:
        cohort[method['name']] = readPredictResult(method['file'])
    ignore_sample = ["HG01123", "HG02486", "HG02559"]
    for s in ignore_sample:
        del cohort["answer"][s]

    result = compare_method(cohort, "answer", "kir", ignore_suffix=True)
    # print(result)
    result_df = result.to_dataframe()


    summary_dfs = []
    for remove_res in [0,3,5,7]:
        result_no_novel_df = result_df.groupby("id", as_index=False).apply(partial(removeNovel, novel_res=remove_res))
        summary = table_summarize(result_no_novel_df)
        summary["remove_res"] = remove_res
        summary_dfs.append(summary)
        with pd.option_context('display.precision', 3):
            print(compact_summary(summary))

    summary_df = pd.concat(summary_dfs)
    with pd.option_context('display.precision', 3):
        print(compact_summary(summary_df, group_by=["remove_res", "method"]))
