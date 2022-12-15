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
    allele_fields = allele_fields.str.findall("(KIR.*?)\(").explode()
    allele_fields = allele_fields[allele_fields.notna()]
    # print(allele_fields)
    return pd.Series([list(allele_fields)], dtype="object")


def mergeHaplo(sample: pd.DataFrame) -> pd.Series:
    # sort by KIR name, but the order will not preserved
    sample["alleles"] = sample["alleles"].apply(sorted)
    return pd.Series({
        "haplos": "+".join(sample["alleles"].str.len().astype(str)),
        "alleles": "_".join(sample["alleles"].sum()),
    })


if __name__ == "__main__":
    hprc_kir_df = pd.read_csv("kelvin_kir.csv", header=None)
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
    print(hprc_kir_df_sample)
    hprc_kir_df_sample.to_csv("hprc_summary.csv", sep="\t", index=False)
