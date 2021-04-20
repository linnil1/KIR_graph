"""
Author: linnil1
Description: Read IPD-KIR to MSA format for furthur development
"""
import re
import os
import copy
import logging
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from pyHLAMSA import KIRmsa, Genemsa

# setup logger to stdout
logger = logging.getLogger("pyHLAMSA")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
logger.addHandler(ch)
thread = 30
# KIR configuration
kir_column = ["5UTR", "exon1", "intron1", "exon2", "intron2", "exon3",
              "intron3", "exon4", "intron4", "exon5", "intron5", "exon6",
              "intron6", "exon7", "intron7", "exon8", "intron8", "exon9",
              "3UTR"]


def run_dk(cmd):
    """ run docker container """
    run("docker run -it --rm --security-opt label=disable -u root -w /app "
        "-v $PWD:/app " + cmd)


def run(cmd):
    """ wrap os.system """
    print(cmd)
    os.system(cmd)


def extract_main_allele(gene):
    """ consensus among each gene sub-group """
    gene_name = gene.gene_name
    gene_sh = Genemsa(gene_name, "gen",
                      gene.blocks, gene.labels)

    consenus_names = map(lambda i: i.split("*")[0] + "*" + i.split("*")[1][:3],
                         gene.alleles.keys())

    for consenus_name in sorted(set(consenus_names)):
        print(consenus_name)
        consenus_allele = gene.select_allele(consenus_name.replace("*", r"\*")
                                             + ".*")
        seq = consenus_allele.get_consensus(include_gap=True)
        gene_sh.add(consenus_name, seq)

    return gene_sh


def align_genes(genes):
    """ Align all KIR genes together """
    genes = copy.deepcopy(genes)
    for gene_name, gene in genes.items():
        label = [("intron", "intron3/4"),  # 2DL4-5
                 ("intron", "intron5/6")]  # 3DL3
        for lab in label:
            if lab not in gene.labels:
                continue
            print(gene_name, "has imcomplete sequence in", lab)
            ind = gene.labels.index(lab)
            i = int(re.findall(r"\d+", lab[1])[0])
            gene.labels[ind] = ("intron", f"intron{i}")
            gene.labels[ind+1:ind+1] = [("exon", f"exon{i+1}"),
                                        ("intron", f"intron{i+1}")]
            gene.blocks[ind+1:ind+1] = [0, 0]

        # 3DP1
        for name in kir_column:
            if name not in map(lambda i: i[1], gene.labels):
                gene.labels.insert(-1, ("intron", name))
                gene.blocks.insert(-1, 0)
                print(gene_name, "has imcomplete sequence in", name)

    return genes


def summary(genes):
    """ Output intron/exon length for each gene """
    df = pd.DataFrame(columns=["name"] + kir_column)
    for gene_name, gene in genes.items():
        gdic = dict(zip(list(map(lambda i: i[1], gene.labels)), gene.blocks))
        gdic['name'] = gene_name
        df = df.append(gdic, ignore_index=True)

    print(df)
    df = df.astype(int, errors="ignore")
    return df


def clustalo(name):
    """ Run clustalomega """
    print("clustalomega alignment")
    run_dk("quay.io/biocontainers/clustalo:1.2.4--h1b792b2_4 "
           f"clustalo --infile {name}.fa -o {name}.aln.fa "
           f"--outfmt fasta --threads {thread}")


def msa_for_chunk(kir_chunk):
    """ Run clustalomega for each intron and exon """
    # write to file
    for chunk_name, chunk in kir_chunk.items():
        name = f"KIR_{chunk_name}"
        with open(name + ".fa", "w") as output_handle:
            chunk_noempty = filter(lambda i: len(i.seq), chunk)
            SeqIO.write(chunk_noempty, output_handle, "fasta")
            # SeqIO.write(chunk, output_handle, "fasta")

    # align
    with ThreadPoolExecutor(max_workers=thread) as executor:
        for chunk_name, chunk in kir_chunk.items():
            name = f"KIR_{chunk_name}"
            executor.submit(clustalo, name)

    # read MSA
    kir_align = {}
    for chunk_name, chunk in kir_chunk.items():
        name = f"KIR_{chunk_name}.aln"
        chunk_new = AlignIO.read(name + ".fa", "fasta")
        length = chunk_new.get_alignment_length()
        chunk_empty = [i.id for i in chunk if not len(i.seq)]
        for empty_chunk_name in chunk_empty:
            chunk_new.append(SeqRecord(seq="-" * length,
                                       name=empty_chunk_name,
                                       id=empty_chunk_name))
        chunk_new.sort()
        kir_align[chunk_name] = chunk_new

    return kir_align


def kir_to_msa():
    """ Main function """
    # read
    kir = KIRmsa(filetype=["gen"])

    # shrink
    kir_sh = {}
    for gene_name in kir.genes:
        kir_sh[gene_name] = extract_main_allele(kir.genes[gene_name])
    summary(kir_sh)

    # align for gene, fill the imcomplete
    kir_draft = align_genes(kir_sh)
    summary(kir_draft)

    # separate chunk
    kir_chunk = {}
    for chunk_i, chunk_name in enumerate(kir_column):
        kir_chunk[chunk_name] = []
        for gene_name in kir_draft:
            kir_chunk[chunk_name].extend(
                    kir_draft[gene_name].select_chunk([chunk_i])
                                        .to_fasta(gap=False))

    # align for each chunk
    kir_align = msa_for_chunk(kir_chunk)

    # merge the chunk
    kir_biomsa = None
    blocks = []
    for chunk_name in kir_column:
        blocks.append(kir_align[chunk_name].get_alignment_length())
        if kir_biomsa is not None:
            kir_biomsa = kir_biomsa + kir_align[chunk_name]
        else:
            kir_biomsa = kir_align[chunk_name]

    # MultipleSeqAlignment -> GeneMSA
    SeqIO.write(kir_biomsa, "kir_merge.fa", "fasta")
    kir_msa_align = Genemsa.from_MultipleSeqAlignment(kir_biomsa)
    kir_msa_align.blocks = blocks
    kir_msa_align.labels = next(iter(kir_draft.values())).labels

    # save to other types
    kir_msa_align.add("KIR*consensus",
                      kir_msa_align.get_consensus(include_gap=False))
    kir_msa_align.save_fasta("kir_merge.consensus.fa", gap=False)
    kir_msa_align.save_bam("kir_merge.bam", "KIR*consensus")
    kir_msa_align.save_gff("kir_merge.gff")
    kir_msa_align.save_msa("kir_merge.save.fa", "kir_merge.save.gff")


def download():
    """ Download needed data """
    run("git clone https://github.com/ANHIG/IPDKIR")
    run("git clone https://github.com/linnil1/pyHLAMSA")
    run("git clone git@github.com:linnil1/pyHLAMSA.git")
    run("git clone git@github.com:linnil1/HLAMSA-export-to-Hisat2.git")
    run("wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads -m -nd -np -P giab_hg002")
    run("cat giab_hg002/D1_S1_L001_R1_0* > merge_D1_S1_L001_R1.fastq.gz")
    run("cat giab_hg002/D1_S1_L001_R2_0* > merge_D1_S1_L001_R2.fastq.gz")


if __name__ == "__main__":
    download()
    kir_to_msa()


"""
Next:
kir_msa = Genemsa.load_msa("kir_merge.save.fa", "kir_merge.save.gff")
python3 kir_to_hisat2.py
/root/hisatgenotype/hisat2/hisat2-build-s --wrapper basic-0 kir_backbone.fa \
        --snp kir.index.snp --haplotype kir.haplotype -p 30 kir.graph
hisat2 -x kir.graph  --no-unal --threads 25 \
    -1 merge_D1_S1_L001_R1.fastq.gz \
    -2 merge_D1_S1_L001_R2.fastq.gz > giab_merge.sam
samtools sort giab_merge.sam -o giab_merge.bam
samtools view giab_merge.bam -f 0x2 -F 256 -o giab_merge.pair.bam
samtools index giab_merge.pair.bam
"""
