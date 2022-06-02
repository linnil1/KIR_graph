import re
import os
import copy
import random
from pathlib import Path
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from pyhlamsa import KIRmsa
from kg_utils import runDocker, runShell


def readHaplo():
    """
    Read the common haplotype configuration

    This csv file is from ping_paper_scripts/synthetic_sequence/KIR_gene_haplotypes.csv
    """
    data = pd.read_csv("KIR_gene_haplotypes.csv")
    data.index = data["hapID"]
    print("Haplotype")
    print(data)
    return data


def readSequences():
    """
    Read sequences

    Return:
      dict[str, str]: sequence-name and sequences
    """
    # any sequences with _sequences.fa are the same
    seqs_dict = SeqIO.to_dict(SeqIO.parse("index/kir_2100_raw.mut01_sequences.fa", "fasta"))
    for i in seqs_dict.values():
        i.description = ""
    return seqs_dict

    kir = KIRmsa(filetype=["gen"], version="2100")
    seqs_dict = {}
    for msa in kir.values():
        seqs_dict.update({
            i.id: i for i in msa.to_fasta(gap=False)
        })
    return seqs_dict


def groupSequences(seqs_dict):
    """
    Group the sequences by gene name

    Args:
      seqs_dict(dict[str, str]): sequence-name and sequences
    Return:
      gene_dict(dict[str, list[Seq]]): gene-name and sequences
    """
    gene_dict = defaultdict(list)
    for i, seq in seqs_dict.items():
        # assume 2DL5A and 2DL5 is same
        gene_dict[i[:7]].append(seq)
    return gene_dict


def randomSelectAlleles(haplo_df, gene_dict):
    """
    Main function, randomly select KIR alleles for sample
    Args:
      haplo_df(DataFrame): HaploType data
      gene_dict(dict[str, Seq]): Key: gene Value: sequences in the gene

    Returns
      haplos(list[str]): Name of haplotypes (Must be 2)
      alleles(list[str]): Name of alleles
    """
    haplos = random.choices(haplo_df.index, k=2)
    gene_count = haplo_df.loc[haplos, :].sum(axis=0)
    gene_count = gene_count.to_dict()
    gene_count = {k: v for k, v in gene_count.items() if k.startswith("KIR")}
    alleles = []
    for gene, count in gene_count.items():
        alleles.extend(random.choices(gene_dict[gene], k=count))
    allele_name_count = defaultdict(int)
    # copy the seq to prevent homozygous edit two time on same object
    alleles = [copy.deepcopy(allele) for allele in alleles]
    for seq in alleles:
        allele_name_count[seq.id] += 1
        seq.id = f"{seq.id}-{allele_name_count[seq.id]}"
    return haplos, alleles


def generateFastq(name, depth=50, fragment_length=400):
    """
    Use art_illumina to generate the fastq

    Run with these parameters
    * -ss HS25  Pair-end
    * -f        depth
    * -m        fragment_length = left + insert + right
    * -na       don't output alignment(.aln)
    * -ef       error free
    * -l        read-length
    * -s        length devation
    * -rs       randomseed
    * -sam      Output sam
    """
    if not os.path.exists("art_illumina"):
        # not exists -> download
        runShell("wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz")
        runShell("tar xf artbinmountrainier2016.06.05linux64.tgz")
        runShell("mv art_bin_MountRainier/art_illumina .")
        runShell("rm art_bin_MountRainier -rf")
        runShell("rm artbinmountrainier2016.06.05linux64.tgz")
    runDocker("alpine", f"""\
        ./art_illumina \
        -ss HS25 \
        -i {name}.fa \
        -l 150 \
        -f {depth} \
        -m {fragment_length} \
        -s 10 \
        -sam -na\
        -rs 444 \
        -o {name}.read. \
    """)


def rewriteSam(file_in, file_out):
    """
    Remove -1 -2 suffix in reference
    """
    s = set()
    with open(file_in) as f_in, open(file_out, "w") as f_out:
        for line in f_in:
            # KIR2DL1*039-1
            line = re.sub(r"KIR(.*?)\*(\w+?)\-\d\t", r"KIR\1*\2\t", line)
            if line in s:
                continue
            s.add(line)
            f_out.write(line)


if __name__ == "__main__":
    N = 10
    # output_name = "linnil1_syn"
    # output_name = "linnil1_syn_full"  # m = 200, N=10, seed=444
    # output_name = "linnil1_syn_wide"  # m = 400, N=100, seed=44
    basename = "linnil1_syn_30x/linnil1_syn_30x"
    depth = 30
    basename = "linnil1_syn_exon/linnil1_syn_exon"
    depth = 90

    # read
    seqs_dict = readSequences()
    seqs_dict = {k: v for k, v in seqs_dict.items() if len(v.seq) > 4200}  # remove exon-only sequence, min 3DP1 4240
    gene_dict = groupSequences(seqs_dict)
    haplo = readHaplo()
    random.seed(44)

    # main
    summary = []
    Path(basename).parent.mkdir(parents=True, exist_ok=True)
    for id in range(N):
        name = f"{basename}.{id:02d}"
        haplos, alleles = randomSelectAlleles(haplo, gene_dict)
        print(name, haplos, [i.id for i in alleles])
        summary.append({
            "id": f"{id:02d}",
            "name": name,
            "haplos": "_".join(haplos),
            "alleles": "_".join([i.id.split("-")[0] for i in alleles]),
        })
        SeqIO.write(alleles, f"{name}.fa", "fasta")
        generateFastq(name, depth=depth)
        rewriteSam(name + ".read..sam", name + ".sam")

    # write summary
    summary = pd.DataFrame(summary)
    print(summary)
    summary.to_csv(f"{basename}.summary.csv", sep="\t", index=False)
