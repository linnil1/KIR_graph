import random
from pathlib import Path

from Bio import SeqRecord, SeqIO
import pandas as pd

from pyhlamsa import KIRmsa, Genemsa

from kg_create_data import readHaplo, randomTwoHaplo

KIR = None


def createRandomFakeIntron(msa: Genemsa) -> tuple[str, str]:
    alleles_full = msa.select_complete().list_alleles()
    alleles_half = msa.select_incomplete().list_alleles()
    allele_full = random.choice(alleles_full)
    allele_half = random.choice(alleles_half)
    msa_part = msa.select_allele([allele_half, allele_full])
    msa_part = msa_part.fill_incomplete(allele_full)
    return f"{allele_half}-intron{allele_full}", msa_part.get(allele_half)


def getRandomFullAllele(msa: Genemsa) -> tuple[str, str]:
    alleles_full = msa.select_complete().list_alleles()
    allele_full = random.choice(alleles_full)
    return allele_full, msa.get(allele_full)


def createFakeIntron(msa: Genemsa, cn: int, fake_num: int = 0) -> Genemsa:
    new_msa = msa.copy(copy_allele=False)
    for i in range(cn):
        if i < fake_num:
            allele, allele_str = createRandomFakeIntron(msa)
        else:
            allele, allele_str = getRandomFullAllele(msa)
        new_msa.append(allele, allele_str)
    return new_msa


def createFakeIntronByCN(
    gene_counts: dict[str, int], fake_num: int = 0
) -> list[SeqRecord]:
    alleles = []
    global KIR
    if not KIR:
        KIR = KIRmsa(filetype=["nuc", "gen"], ipd_folder="KIR_v2100", version="2100")

    for gene, msa in KIR.items():
        if gene_counts.get(gene):
            new_msa = createFakeIntron(msa, gene_counts[gene], fake_num)
            alleles.extend(new_msa.to_records(gap=False))
    return alleles


def createFakeIntronSample(N: int, basename: str, seed: int = 44, fake_num: int = 1):
    """Random pick the alleles for N samples"""
    # read
    haplo_df = readHaplo()
    random.seed(seed)

    # main
    summary = []
    Path(basename).parent.mkdir(parents=True, exist_ok=True)
    for id in range(N):
        name = f"{basename}.{id:02d}"
        haplos, gene_counts = randomTwoHaplo(haplo_df)
        alleles = createFakeIntronByCN(gene_counts, fake_num)
        print(name, haplos, [i.id for i in alleles])
        summary.append(
            {
                "id": f"{id:02d}",
                "name": name,
                "haplos": "_".join(haplos),
                "alleles": "_".join([i.id for i in alleles]),
            }
        )
        SeqIO.write(alleles, f"{name}.fa", "fasta")

    # write summary
    summary = pd.DataFrame(summary)
    print(summary)
    summary.to_csv(f"{basename}_summary.csv", sep="\t", index=False)


if __name__ == "__main__":
    N = 1
    basename = "test/linnil1_syn_fakeintron"
    createFakeIntronSample(N, basename, seed=444)
