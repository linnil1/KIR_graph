import glob
from pathlib import Path
from Bio import SeqIO
from kg_utils import runDocker, runShell, samtobam
from pyhlamsa import KIRmsa, msaio

kir_default = None


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


def rename(name_ori, output_name):
    runShell(f"ln -fs ../{name_ori}.1.fq {output_name}.read.1.fq")
    runShell(f"ln -fs ../{name_ori}.2.fq {output_name}.read.2.fq")
    return output_name


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


def mergeGff(gff_files, name):
    with open(f"{name}.gff", "w") as f_gff:
        f_gff.write("##gff-version 3\n")
        for f_gff_allele in gff_files:
            allele = f_gff_allele.split('.')[-2]
            with open(f_gff_allele) as f:
                for i in f:
                    if i.startswith("##"):
                        continue
                    # rename ref to ref-1 ref-2
                    f_gff.write("\t".join([allele, *i.strip().split("\t")[1:]]) + "\n")
            runShell(f"rm '{f_gff_allele}'")


def extractExonPairReads(reference_fasta, bam_file, output_name, kir=None):
    seqs = SeqIO.parse(reference_fasta, "fasta")
    name = ".".join(bam_file.split('.')[:-1])

    # save kir object in global
    global kir_default
    if not kir:
        if not kir_default:
            kir_default = KIRmsa(filetype=["gen"], version="2100")
        kir = kir_default

    regions = []
    for seq in seqs:
        ref = seq.id[:7]  # KIR2DL5AB -> KIR2DL5
        regions.extend(getExonRegion(kir[ref], seq, name))
    mergeGff(list(glob.glob(f"{name}.KIR*.gff")), name)

    print(regions)
    samtobam(name, keep=True)
    name_overlap = extractRegion(name, regions)
    name_extract = reservePair(name_overlap, name)
    runShell(f"ln -fs ../{name_extract}.bam     {output_name}.bam")
    runShell(f"ln -fs ../{name_extract}.bam.bai {output_name}.bam.bai")
    name_extract = bam2fastq(name_extract)
    name = rename(name_extract, output_name)
    return name


if __name__ == "__main__":
    """
    Run this after kg_create_data.py
    """
    kir = KIRmsa(filetype=["gen"], version="2100")
    for i in range(10):
        name = f"linnil1_exon_30x/linnil1_exon_30x.{i:02d}"
        name = extractExonPairReads(name + ".fa", name + ".read..sam", name + "_exon", kir=kir)
        print(name)
