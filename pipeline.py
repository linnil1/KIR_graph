"""
Author: linnil1
Description: Read IPD-KIR to MSA format for furthur development
"""
import re
import os
import glob
import copy
import logging
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from pyHLAMSA import KIRmsa, Genemsa
from collections import defaultdict

# this require hisat2 docker
import sys
sys.path.append('./hisat-genotype/hisatgenotype_modules')
os.environ["PATH"] = os.environ["PATH"] + ":./"

try:  # temp
    import pipeline_to_hisat
    import pipeline_typing
except:
    pass

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
genes = ["KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5",
         "KIR2DP1", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5",
         "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DP1", "KIR3DS1"]
fastq_suffix = ".fastq"


def run_dk(image, cmd):
    """ run docker container """
    run("docker run -it --rm --security-opt label=disable -u root -w /app "
        f"-v $PWD:/app {image} {cmd}")


def run(cmd):
    """ wrap os.system """
    print(cmd)
    os.system(cmd)


def extract_main_allele(gene, allele_group=True):
    """ consensus among each gene sub-group """
    gene_name = gene.gene_name
    gene_sh = Genemsa(gene_name, "gen",
                      gene.blocks, gene.labels)

    if allele_group:
        consenus_names = map(lambda i: i.split("*")[0] + "*" + i.split("*")[1][:3],
                             gene.alleles.keys())
    else:
        consenus_names = gene.alleles.keys()

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
            """
            # method 1:
            gene.labels[ind] = ("intron", f"intron{i}")
            gene.labels[ind+1:ind+1] = [("exon", f"exon{i+1}"),
                                        ("intron", f"intron{i+1}")]
            gene.blocks[ind+1:ind+1] = [0, 0]
            """
            # method 2:
            gene.labels[ind] = ("intron", f"intron{i+1}")
            gene.labels[ind:ind] = [("intron", f"intron{i}"),
                                    ("exon", f"exon{i+1}")]
            gene.blocks[ind:ind] = [0, 0]

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
    run_dk("quay.io/biocontainers/clustalo:1.2.4--h1b792b2_4",
           f"clustalo --infile {name}.fa -o {name}.aln.fa "
           f"--outfmt fasta --threads {thread} --force")


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
    with ProcessPoolExecutor(max_workers=thread) as executor:
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


def kir_to_msa(name="kir_merge", allele_group=True):
    """ Main function """
    # read
    kir = KIRmsa(filetype=["gen"])

    # shrink
    kir_sh = {}
    for gene_name in kir.genes:
        kir_sh[gene_name] = extract_main_allele(kir.genes[gene_name], allele_group)
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
    SeqIO.write(kir_biomsa, f"{name}.fa", "fasta")
    kir_msa_align = Genemsa.from_MultipleSeqAlignment(kir_biomsa)
    kir_msa_align.blocks = blocks
    kir_msa_align.labels = next(iter(kir_draft.values())).labels

    # save to other types
    kir_msa_align.add("KIR*BACKBONE",
                      kir_msa_align.get_consensus(include_gap=False))
    kir_msa_align.save_fasta(f"{name}.consensus.fa", gap=False)
    kir_msa_align.save_bam(f"{name}.bam", "KIR*BACKBONE")
    kir_msa_align.save_gff(f"{name}.gff")
    kir_msa_align.save_msa(f"{name}.save.fa", f"{name}.save.gff")


def kir_to_multi_msa():
    """ Main function """
    # read
    kir = KIRmsa(filetype=["gen"])

    # shrink
    for gene_name, kir_msa_align in kir.genes.items():
        kir_msa_align.add(f"{gene_name}*BACKBONE",
                          kir_msa_align.get_consensus(include_gap=False))
        kir_msa_align.save_bam(f"kir_split.{gene_name}.bam", f"{gene_name}*BACKBONE")
        kir_msa_align.save_gff(f"kir_split.{gene_name}.gff")
        kir_msa_align.save_msa(f"kir_split.{gene_name}.save.fa", f"kir_split.{gene_name}.save.gff")


def download():
    """ Download needed data """
    run("git clone https://github.com/ANHIG/IPDKIR")
    run("git clone https://github.com/linnil1/pyHLAMSA")
    run("git clone git@github.com:linnil1/pyHLAMSA.git")
    run("git clone git@github.com:linnil1/HLAMSA-export-to-Hisat2.git")


def download_data():
    """ Download large data """
    run("wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads -m -nd -np -P giab_hg002")
    run("cat giab_hg002/D1_S1_L001_R1_0* > merge_D1_S1_L001_R1.fastq.gz")
    run("cat giab_hg002/D1_S1_L001_R2_0* > merge_D1_S1_L001_R2.fastq.gz")



def hisat2(index="kir_merge"):
    for name in samples:
        f1, f2 = name + ".read1" + fastq_suffix, name + ".read2" + fastq_suffix
        # run(f"hisat2 --threads 25 -x {index}.graph -1 {f1} -2 {f2} --no-unal --no-spliced-alignment --max-altstried 64 --haplotype > {name}.sam")
        run(f""" \
            hisat2 --threads {thread} -x {index}.graph -1 {f1} -2 {f2} \
            --no-spliced-alignment --max-altstried 64 --haplotype \
            > {name}{suffix}.sam
        """)


def bowtie2(index="kir_split.linear"):
    # dk quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0
    for name in samples:
        f1, f2 = name + ".read1" + fastq_suffix, name + ".read2" + fastq_suffix
        # run(f"bowtie2 --threads 25 -x {index}.linear -1 {f1} -2 {f2} --no-unal  > {name}.linear.sam")
        run_dk("quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0",
               f"bowtie2 --threads 25 -x {index} -1 {f1} -2 {f2} -a -S {name}{suffix}.sam")


def bowtie2Ping():
    index = "kir_split.linear"
    for name in samples:
        f1, f2 = name + ".read1" + fastq_suffix, name + ".read2" + fastq_suffix
        # run(f"bowtie2 --threads 25 -x {index}.linear -1 {f1} -2 {f2} --no-unal  > {name}.linear.sam")
        # '--no-unal',
        run(f"bowtie2 --threads 25 -x {index} -1 {f1} -2 {f2} " + \
            " ".join(['-5 0', '-3 6', '-N 0', '--end-to-end', '--score-min L,-2,-0.08', '-I 75', '-X 1000', '-a','--np 1', '--mp 2,2', '--rdg 1,1', '--rfg 1,1']) + \
            f"-S  {name}{suffix}.sam")


def bowtie2BuildFull():
    index = "kir_full"
    run_dk("quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0",
           f"bowtie2-build kir_merge_sequences.fa {index} --threads 30")


def bowtie2Full():
    index = "kir_full"
    for name in samples:
        f1, f2 = name + ".read1" + fastq_suffix, name + ".read2" + fastq_suffix
        # '--no-unal',
        run_dk("quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0",
               f"bowtie2 --threads 25 -x {index} -1 {f1} -2 {f2} -X 1000 -a --end-to-end -S {name}.full.sam")


def bowtie2Build():
    index = "kir_split.linear"
    run_dk("quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0",
           f"bowtie2-build kir_split_backbone.fa {index} --threads 30")

def bowtie2BuildNotgroup():
    index = "kir_split_full.linear"
    run_dk("quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0",
           f"bowtie2-build kir_split_full_backbone.fa {index} --threads 30")

def bowtie2BuildFull():
    index = "kir_split_full.linear"
    run_dk("quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0",
           f"bowtie2-build kir_split_full_backbone.fa {index} --threads 30")


def samtobam():
    for name in samples:
        name += suffix
        run(f"samtools sort {name}.sam -o {name}.bam")
        run(f"samtools view {name}.bam -f 0x2 -F 256 -o {name}.pair.bam")
        run(f"samtools index {name}.pair.bam")


def fastqc():
    run_dk("docker.io/biocontainers/fastqc:v0.11.9_cv7",
           f"clustalo --infile {name}.fa -o {name}.aln.fa "
           f"--outfmt fasta --threads {thread} --force")


def hisatTyping(index):
    pipeline_typing.hisatdataInit(index)
    for sample in samples:
        pipeline_typing.typingAllGene(f"{sample}{suffix}.pair.bam")
    # with ProcessPoolExecutor(max_workers=thread) as executor:
    #     for sample in samples:
    #         executor.submit(pipeline_typing.typingAllGene, f"{sample}{suffix}.pair.bam")


def linkFastq():
    for i in range(10):
        name = "linnil1_syn_full"
        os.system(f"ln -fs ../{name}/{name}.{i:02d}.read.1.fq data/{name}.{i:02d}.read1.fastq")
        os.system(f"ln -fs ../{name}/{name}.{i:02d}.read.2.fq data/{name}.{i:02d}.read2.fastq")


if __name__ == "__main__":
    # download()
    # download_data()
    samples = sorted(map(lambda i: i.split(".read")[0], glob.glob("data/synSeq.*.read1.fastq")))

    # create my sample
    # python3 pipeline_generate_syn.py
    samples = ["data/synSeq.hiseq.dp50.rl150.1"]
    samples = [f"data/linnil1_syn.0{i}" for i in range(10)]
    samples = [f"data/linnil1_syn_full.0{i}" for i in range(10)]
    samples = samples[:1]

    # kir_merge
    # kir_to_msa()
    # pipeline_to_hisat.main("kir_merge")
    # pipeline_to_hisat.build("kir_merge")
    suffix = ".merge"
    # hisat2()
    # samtobam()
    # hisatTyping("./kir_merge")

    # kir_split
    # kir_to_multi_msa()
    # pipeline_to_hisat.main("kir_split")
    # pipeline_to_hisat.build("kir_split")
    # run("cat kir_split.*.save.gff > kir_split.gff")

    suffix = ".linear"
    # bowtie2Build()
    # bowtie2()
    # samtobam()

    suffix = ".full"
    # bowtie2BuildFull()
    # bowtie2Full()
    # samtobam()

    suffix = ".split"
    # hisat2("kir_split")
    # samtobam()
    # pipeline_typing.hisatdataInit("./kir_split")
    # pipeline_typing.typingAllGene(f"{samples[0]}{suffix}.pair.bam")

    suffix = ".merge"
    # kir_to_msa("kir_merge_full", allele_group=False)
    # pipeline_to_hisat.main("kir_merge_full")
    # pipeline_to_hisat.build("kir_merge_full")
    # os.system("python3 pipeline_generate_syn.py")
    # hisat2("kir_merge_full")
    # samtobam()
    hisatTyping("./kir_merge_full")

    suffix = ".split"
    # pipeline_to_hisat.main("kir_split_full")
    # pipeline_to_hisat.build("kir_split_full")
    # hisat2("kir_split_full")
    # samtobam()
    hisatTyping("./kir_split_full")

    suffix = ".linear"
    # bowtie2BuildNotgroup()
    # bowtie2()
    # samtobam()

"""
Tips:
# load data(In script)
kir_msa = Genemsa.load_msa("kir_merge.save.fa", "kir_merge.save.gff")

# MSA to hisat index
python3 kir_to_hisat2.py
/root/hisatgenotype/hisat2/hisat2-build-s --wrapper basic-0 kir_backbone.fa \
        --snp kir.index.snp --haplotype kir.haplotype -p 30 kir.graph

# Read Mapping: Sometimes Segment Error(I don't know why)
hisat2 -x kir.graph --threads 25 --no-unal --no-spliced-alignment \
        --max-altstried 64 --haplotype -X 2000 \
        -1 merge_D1_S1_L001_R1.fastq.gz \
        -2 merge_D1_S1_L001_R2.fastq.gz > giab_merge.sam
samtools sort giab_merge.sam -o giab_merge.bam

# Remove singleton
samtools view giab_merge.bam -f 0x2 -F 256 -o giab_merge.pair.bam
samtools index giab_merge.pair.bam

# Remove LINE region
echo KIR*consensus$'\t'18400$'\t'19700 > kir.bed
bedtools intersect -a giab_merge.pair.bam -b ../kir.bed -v | samtools sort - -o giab_merge.pair.filter.bam
samtools index giab_merge.filter.bam
"""
