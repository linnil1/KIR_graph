from pathlib import Path
from collections import Counter

from graphkir.utils import (
    samtobam,
    getThreads,
)
from graphkir.hisat2 import (
    hisatMap,
    readPair,
    filterRead,
    ReadsAndVariantsData,
    PairRead,
    saveReadsToBam,
)
from graphkir import wgs

from kg_utils import runDocker


def hisatMapWrap(input_name, index):
    # 1 to 1
    output_name = input_name + "." + index.replace("/", "_")
    if Path(f"{output_name}.bam").exists():
        return output_name
    f1, f2 = input_name + ".read.1.fq.gz", input_name + ".read.2.fq.gz"
    if not Path(f1).exists():
        f1, f2 = input_name + ".read.1.fq", input_name + ".read.2.fq"
    hisatMap(index, f1, f2, output_name + ".bam", threads=getThreads())
    return output_name


def bwaIndex(input_name="index/kir_2100_raw.mut01"):
    # brute force (remove mut01)
    output_name = input_name + ".bwa"
    if Path(output_name + ".bwt").exists():
        return output_name
    if Path(f"{input_name}.fa").exists():
        fasta = f"{input_name}.fa"
    elif Path(f"{input_name}.fa.gz").exists():
        fasta = f"{input_name}.fa.gz"
    else:
        raise ValueError("fasta not found: " + str(input_name))
    wgs.bwaIndex(fasta, output_name)
    return output_name


def bwa(input_name, index):
    suffix = "." + index.replace("/", "_")
    args = ""
    output_name = input_name + suffix
    if Path(f"{output_name}.bam").exists():
        return output_name

    # main
    f1, f2 = input_name + ".read.1.fq", input_name + ".read.2.fq"
    if not Path(f1).exists():
        f1, f2 = input_name + ".read.1.fq.gz", input_name + ".read.2.fq.gz"
    if not Path(f1).exists():
        raise ValueError("fastq not found: " + str(input_name) + ".read.{}.fq")
    wgs.bwa(index, f1, f2, output_name, threads=getThreads())
    return output_name


def bowtie2Index(input_name):
    output_name = input_name + ".bowtie2"
    if Path(output_name + ".1.bt2").exists():
        return output_name
    runDocker("bowtie",
              f"bowtie2-build {input_name}.fa {output_name} --threads {getThreads()}")
    return output_name


def bowtie2(input_name, index, use_arg="default"):
    suffix = "." + index.replace("/", "_")
    args = ""
    if use_arg == "ping":
        args = "  -5 0  -3 6  -N 0  --end-to-end  --score-min L,-2,-0.08  " + \
               "  -I 75  -X 1000  -a  --np 1  --mp 2,2  --rdg 1,1  --rfg 1,1  "
        suffix += "_ping"

    output_name = input_name + suffix
    if Path(f"{output_name}.bam").exists():
        return output_name

    # main
    f1, f2 = input_name + ".read.1.fq", input_name + ".read.2.fq"
    runDocker("bowtie",
              f"bowtie2 {args} --threads {getThreads()} -x {index} "
              f"-1 {f1} -2 {f2} -a -S {output_name}.sam")
    samtobam(output_name)
    return output_name


def trimBam(input_name, remove=False):
    """ Trim very large bam """
    output_name = input_name + ".trim"
    if Path(f"{output_name}.bam").exists():
        return output_name
    runDocker("samtools", f"samtools view -h -@ {getThreads()} -F 4 {input_name}.bam -o {output_name}.sam")
    samtobam(output_name)
    if remove:
        Path(f"{input_name}.bam").unlink()
        Path(f"{input_name}.bam").touch()
    return output_name


def filterNM(input_name):
    # name = "data/linnil1_syn_s2022.99.30x_s1031.index_kir_2100_withexon_ab_2dl1s1.leftalign.backbone.bwa"
    output_prefix = input_name + ".NM_4"
    bam_file = input_name + ".bam"
    if Path(output_prefix + ".no_multi.bam").exists():
        return output_prefix + ".no_multi"
    pair_reads = readPair(bam_file)
    pair_reads = filter(lambda lr: filterRead(lr[0]) and filterRead(lr[1]), pair_reads)
    pair_reads = list(pair_reads)

    c = Counter(left.split("\t")[0] for left, right in pair_reads)
    for i, j in c.items():
        if j > 1:
            print(i, j)

    reads = []
    for left_record, right_record in pair_reads:
        reads.append(
            PairRead(
                lpv=[],
                lnv=[],
                rpv=[],
                rnv=[],
                l_sam=left_record,
                r_sam=right_record,
                multiple=c[left_record.split("\t")[0]],
                backbone=left_record.split("\t")[2],
            )
        )

    reads_data: ReadsAndVariantsData = {
        "variants": [],
        "reads": reads,
    }
    saveReadsToBam(
        reads_data, output_prefix + ".no_multi", bam_file, filter_multi_mapped=True
    )
    return output_prefix + ".no_multi"
