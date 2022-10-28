import glob
from pathlib import Path
from Bio import SeqIO

from pyhlamsa import KIRmsa, msaio
from graphkir.utils import runShell, samtobam
from kg_utils import runDocker


kir_default = None


def extractRegion(input_name, bed_file):
    output_name = input_name + ".extract"
    if Path(output_name + ".bam").exists():
        return output_name
    runDocker("samtools", f"samtools view -L {bed_file} {input_name}.bam -o {output_name}.bam")
    return output_name


def reservePair(name_filtered, name_original):
    # retrieve all filter-PASSed read id
    proc = runDocker("samtools", f"samtools view {name_filtered}.bam", capture_output=True)
    reads = set([line.split("\t")[0] for line in proc.stdout.split("\n")])

    # filter read if it is in the id set
    proc = runDocker("samtools", f"samtools view {name_original}.bam -h", capture_output=True)
    first_line = True
    with open(f"{name_filtered}.reserve_pair.sam", "w") as f:
        for line in proc.stdout.split("\n"):
            if line.startswith("@") or line.split("\t")[0] in reads:
                if not first_line:
                    f.write("\n")
                f.write(line)
                first_line = False
    name_filtered += ".reserve_pair"
    samtobam(name_filtered)
    return name_filtered


def bam2fastq(name):
    if Path(f"{name}.sortn.read.1.fq").exists():
        return name + ".sortn.read"
    runDocker("samtools", f"samtools sort  -n {name}.bam -o {name}.sortn.bam")
    name += ".sortn"
    runDocker("samtools", f"samtools fastq -n {name}.bam -1 {name}.read.1.fq -2 {name}.read.2.fq -0 /dev/null -s /dev/null")
    name += ".read"
    return name


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
        if pos[0] == pos[1]:
            continue
        yield seq.id, pos[0] + 1, pos[1]


def calcExonToBed(reference_name, kir=None):
    # save kir object in global
    output_name = reference_name + ".exon_region"
    if Path(output_name + ".bed").exists():
        return output_name

    global kir_default
    if not kir:
        if not kir_default:
            kir_default = KIRmsa(filetype=["gen"], version="2100")
        kir = kir_default

    seqs = SeqIO.parse(reference_name + ".fa", "fasta")
    regions = []
    for seq in seqs:
        ref = seq.id[:7]  # KIR2DL5AB -> KIR2DL5
        regions.extend(getExonRegion(kir[ref], seq))

    with open(f"{output_name}.bed", "w") as f:
        for chrom, start, end in regions:
            f.write(f"{chrom}\t{start}\t{end}\n")
    return output_name


def extractPairReadsOnceInBed(input_name, bed_name, kir=None):
    output_name = input_name + ".extract.reserve_pair"
    bed_file = bed_name.format(input_name.template_args[0]) + ".bed"
    if Path(output_name + ".bam").exists():
        return output_name
    bam_bed = extractRegion(input_name, bed_file)
    bam_pair = reservePair(bam_bed, input_name)
    assert str(bam_pair) == str(output_name)
    return output_name
