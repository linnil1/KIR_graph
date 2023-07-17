"""
HPRC's KIR alleles format called from long read by Kelvin.
I transform it to our format
"""
from typing import Any
import pandas as pd


def mergeAlleleField(sample: pd.Series) -> pd.Series:
    # print(sample)
    allele_fields = sample.iloc[3:]
    # empty -> skip
    allele_fields = allele_fields[allele_fields.notna()]
    # KIR3DL3*00101(100.00%) -> KIR3DL3*0010101
    # KIR2DL4*01201(99.91%)~KIR3DL1*03501(100.00%)~KIR2DS4*00103(96.10%)~
    # In version 1.0
    # many characters are introduced
    # KIR2DS4*00101_e124567-3DL1*03501_e8
    # KIR3DL3# (99.92% of 00902)
    # KIR3DL3*01502$
    # KIR3DL3*027=
    # KIR3DL3*035+
    # previous regex
    # allele_fields = allele_fields.str.findall(r"(KIR.*?)\(").explode()
    allele_fields = allele_fields.str.findall(r"KIR[^~^(^+^#^$^=]*").explode()
    allele_fields = allele_fields[allele_fields.notna()]
    allele_fields = allele_fields.str.replace("_", "")
    # print(allele_fields)
    return pd.Series([list(allele_fields)], dtype="object")


def ambiguousTo1(alleles: list[str]) -> list[str]:
    """
    Handle the ambiguous alleles
    origin: KIR2DP1*010/KIR2DP1*006(99.81%)
    after:  KIR2DP1*010
    """
    return [allele.split("/")[0] for allele in alleles]


def mergeHaplo(sample: pd.DataFrame) -> pd.Series:
    # sort by KIR name, but the order will not preserved
    sample["alleles"] = sample["alleles"].apply(sorted)
    return pd.Series({
        "haplos": "+".join(sample["alleles"].str.len().astype(str)),
        "alleles": "_".join(ambiguousTo1(sample["alleles"].sum())),
    })


if __name__ == "__main__":
    input_csv = "kelvin_kir.csv"
    output_csv = "hprc_summary.tsv"
    input_csv = "kelvin_kir_v0_2.csv"
    output_csv = "hprc_summary_v0_2.tsv"
    input_csv = "kelvin_kir_v1_0.csv"
    output_csv = "hprc_summary_v1_0.tsv"
    hprc_kir_df = pd.read_csv(input_csv, header=None)
    if "v1_0" in input_csv:
        hprc_kir_df.loc[1:6, 0] = "hg002"  # version1.0 HG002.hifiasm.hic.0.16.1.rep4.hap1 = HG002
        hprc_kir_df = hprc_kir_df.drop(index=[2,3,4,5])  # hg002 extra

    hprc_kir_df["alleles"] = hprc_kir_df.apply(mergeAlleleField, axis=1)
    print(hprc_kir_df)
    hprc_kir_df = hprc_kir_df[hprc_kir_df["alleles"].str.len() > 0]
    hprc_kir_df = hprc_kir_df.rename(columns={0: "id", 1: "PM"})
    hprc_kir_df = hprc_kir_df[["id", "PM", "alleles"]]
    # hprc_kir_df = hprc_kir_df[hprc_kir_df["id"] != "HG002"]
    print(hprc_kir_df)

    # transform to our format
    hprc_kir_df_sample = hprc_kir_df.groupby("id", as_index=False).apply(mergeHaplo)
    hprc_kir_df_sample["name"] = hprc_kir_df_sample["id"]
    # version 0.2 why ID become lower case
    hprc_kir_df_sample["id"] = hprc_kir_df_sample["id"].str.upper()
    hprc_kir_df_sample["name"] = hprc_kir_df_sample["name"].str.upper()
    print(hprc_kir_df_sample)
    hprc_kir_df_sample.to_csv(output_csv, sep="\t", index=False)
