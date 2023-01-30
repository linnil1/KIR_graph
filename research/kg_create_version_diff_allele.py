from typing import Iterable
from pathlib import Path
from collections import defaultdict
import random

from Bio import SeqIO, SeqRecord
from pyhlamsa import Genemsa
import pandas as pd

from graphkir.utils import readFromMSAs


def getExonOnly(allele_names: Iterable[str]) -> set[str]:
    return set(filter(lambda i: i.endswith("e"), allele_names))


def getFullOnly(allele_names: Iterable[str]) -> set[str]:
    return set(filter(lambda i: not i.endswith("e"), allele_names))


def calcStat(kir_old: dict[str, Genemsa], kir_new: dict[str, Genemsa]) -> None:
    stat = []
    for gene in kir_old:
        allele_del = kir_old[gene].list_alleles() - kir_new[gene].list_alleles()
        print("old-exon -> new-full", allele_del)
        # allele_new = kir_new[gene].list_alleles() - kir_old[gene].list_alleles()
        stat.append(
            {
                "gene": gene,
                "old-exon": len(getExonOnly(kir_old[gene].list_alleles())),
                "old-full": len(getFullOnly(kir_old[gene].list_alleles())),
                "new-exon": len(getExonOnly(kir_new[gene].list_alleles())),
                "new-full": len(getFullOnly(kir_new[gene].list_alleles())),
                "old-exon->new-full": len(getExonOnly(allele_del)),
            }
        )
    df = pd.DataFrame(stat)
    df = df.sort_values("gene")
    df["old-total"] = df["old-exon"] + df["old-full"]
    df["new-total"] = df["new-exon"] + df["new-full"]
    df["exon+"]     = df["new-exon"] - df["old-exon"]
    df["full+"]     = df["new-full"] - df["old-full"]
    df_sum = df.sum(axis=0)
    df_sum["gene"] = "Total"
    df = df.append(df_sum, ignore_index=True)
    print(df)


def randomSelectOldNewAllele(
    kir_old: dict[str, Genemsa],
    kir_new: dict[str, Genemsa],
    allele_use_new: list[bool] = [False, True],
) -> list[SeqRecord]:
    new_alleles = []
    for gene in kir_new:
        old_allele = getFullOnly(kir_old[gene].list_alleles())
        new_allele = getFullOnly(kir_new[gene].list_alleles()
                               - kir_old[gene].list_alleles())
        if not len(new_allele):
            if gene == "KIR3DL3":
                raise ValueError("KIR3DL3 doesn't have new allele")
            continue
        for use_new in allele_use_new:
            if use_new:
                a = random.choices(list(new_allele), k=1)
                msa_a = kir_new[gene].select_allele(a)
            else:
                a = random.choices(list(old_allele), k=1)
                msa_a = kir_old[gene].select_allele(a)
            new_alleles.extend(msa_a.to_records(gap=False))

    # add suffix to prevent duplication
    name_set: dict[str, int] = defaultdict(lambda: 1)
    for allele in new_alleles:
        allele.id = allele.id + "-" + str(name_set[allele.id])
        name_set[allele.id] += 1
    return new_alleles


def removeStrangeAllele(msas: dict[str, Genemsa]) -> None:
    for gene, msa in msas.items():
        allele_del = [gene + "*BACKBONE"]
        for allele, seq in msa.items():
            if len(seq.replace("-", "")) < 3000:
                print("Remove", allele, len(seq.replace("-", "")))
                allele_del.append(allele)
        msa.remove_allele(allele_del)


def createOldNewAllele(input_name: str, N: int, new_allele: int = 1, seed: int = 44) -> str:
    basename = input_name + "/" + input_name + "_old290new2100"
    if new_allele != 1:
        basename += f"_old{new_allele}"
    basename += f"_s{seed}"
    output_name = basename + ".{}"
    if Path(basename + "_summary.csv").exists():
        return output_name

    kir_old_version = readFromMSAs("index/kir_290_withexon")
    kir_new_version = readFromMSAs("index/kir_2100_withexon")
    # kir_old_version = readFromMSAs("index/kir_290_withexon")
    # kir_new_version = readFromMSAs("index/kir_2120_withexon")
    calcStat(kir_old_version, kir_new_version)

    # remove alleles
    removeStrangeAllele(kir_old_version)
    removeStrangeAllele(kir_new_version)

    summary = []
    random.seed(seed)
    for id in range(N):
        name = f"{basename}.{id:02d}"
        if new_allele == 1:
            alleles = randomSelectOldNewAllele(kir_old_version, kir_new_version)
        elif new_allele == 0:
            alleles = randomSelectOldNewAllele(
                kir_old_version, kir_new_version, [False, False]
            )
        else:
            raise ValueError
        print(i.id for i in alleles)
        summary.append(
            {
                "id": f"{id:02d}",
                "name": name,
                "haplos": "1old_1new",
                "alleles": "_".join([i.id.split("-")[0] for i in alleles]),
            }
        )
        SeqIO.write(alleles, f"{name}.fa", "fasta")

    # write summary
    df_summary = pd.DataFrame(summary)
    # print(summary)
    df_summary.to_csv(f"{basename}_summary.csv", sep="\t", index=False)
    return output_name


if __name__ == "__main__":
    N = 1
    createOldNewAllele("test/linnil1_syn_oldnew", N, 1)
