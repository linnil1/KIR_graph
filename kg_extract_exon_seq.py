import os
import glob
from Bio import SeqIO
from kg_utils import runDocker, runShell, samtobam
from pyhlamsa import KIRmsa, msaio


def extractRegion(name, regions):
    runDocker("samtools", f"samtools view {name}.bam -o {name}.exon.bam {' '.join(regions)}")
    return name + ".exon"


def reservePair(name_filtered, name_original):
    proc = runDocker("samtools", f"samtools view {name_filtered}.bam", capture_output=True)
    reads = set([line.split("\t")[0] for line in proc.stdout.split("\n")])
    print(reads)

    proc = runDocker("samtools", f"samtools view {name_original}.bam -h", capture_output=True)
    first = True
    with open(f"{name_filtered}.reserve_pair.sam", "w") as f:
        for line in proc.stdout.split("\n"):
            if line.startswith("@") or line.split("\t")[0] in reads:
                if not first:
                    f.write("\n")
                f.write(line)
                first = False
    name_filtered += ".reserve_pair"
    samtobam(name_filtered)
    return name_filtered


def bam2fastq(name):
    runDocker("samtools", f"samtools sort -n {name}.bam -o {name}.sortn.bam")
    name += ".sortn"
    runDocker("samtools", f"samtools fastq -n {name}.bam -1 {name}.read.1.fq -2 {name}.read.2.fq -0 /dev/null -s /dev/null")
    name += ".read"
    return name


def rename(name, name_ori):
    runShell(f"ln -fs ../{name_ori}.1.fq {name}.exon.read.1.fq")
    runShell(f"ln -fs ../{name_ori}.2.fq {name}.exon.read.2.fq")
    return name + ".exon"


def getExonRegion(msa, seq, gff_name=""):
    allele = seq.id.split("-")[0]
    submsa = msa.select_allele([allele]).shrink().reset_index()  # we need gapless position
    if gff_name:
        msaio.to_gff(submsa, f"{gff_name}.{seq.id}.gff")
    assert len(submsa.alleles[allele]) == len(seq.seq)

    current_pos = 0
    exon_pos = []
    for b in submsa.blocks:
        if b.type == "exon":
            exon_pos.append((current_pos, current_pos + b.length))
        current_pos += b.length
    assert len(submsa.alleles[allele]) == current_pos

    for pos in exon_pos:
        # bed-format = included start + included end
        yield f"{seq.id}:{pos[0]+1}-{pos[1]}"


def mergeGff(name):
    with open(f"{name}.gff", "w") as f_gff:
        f_gff.write("##gff-version 3\n")
        for f_gff_allele in glob.glob(f"{name}.KIR*.gff"):
            allele = f_gff_allele.split('.')[-2]
            with open(f_gff_allele) as f:
                for i in f:
                    if i.startswith("##"):
                        continue
                    # rename ref to ref-1 ref-2
                    f_gff.write("\t".join([allele, *i.strip().split("\t")[1:]]) + "\n")
            runShell(f"rm {f_gff_allele}")


def extractExonPairReads(name, kir=None):
    name_read = name + ".read."
    seqs = SeqIO.parse(f"{name}.fa", "fasta")

    if not kir:
        kir = KIRmsa(filetype=["gen"], version="2100")
    regions = []
    for seq in seqs:
        ref = seq.id[:7]  # AB is same
        regions.extend(getExonRegion(kir[ref], seq, name))
    mergeGff(name)

    print(regions)
    samtobam(name_read, keep=True)
    name_overlap = extractRegion(name_read, regions)
    name_extract = reservePair(name_overlap, name_read)
    name_extract = bam2fastq(name_extract)
    name = rename(name, name_extract)
    return name


if __name__ == "__main__":
    """
    Run this after kg_create_data.py
    """
    kir = KIRmsa(filetype=["gen"], version="2100")
    for i in range(10):
        name = f"linnil1_syn_exon/linnil1_syn_exon.{i:02d}"
        name = extractExonPairReads(name, kir=kir)
        print(name)
