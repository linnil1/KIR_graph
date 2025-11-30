"""
WGS index/mapping part of graphkir
"""

import json

from .utils import getThreads, runShell, samtobam, logger, downloadFile
from .external_tools import runTool
from .samtools_utils import bam2Depth, readSamtoolsDepth

# hs37d5's diploid regions
regions_of_diploid = {
    "VDR": "12:48235320-48298777",
    "RYR1": "19:38924331-39078204",
    "EGFR": "7:55086710-55279321",
}

def downloadHg19(index_folder: str) -> str:
    """Download hs37d5"""
    output_name = f"{index_folder}/hs37d5.fa.gz"
    logger.info(f"[WGS] Download {output_name}")
    url = "https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
    downloadFile(url, output_name)
    return output_name


def bwaIndex(fasta: str, output_name: str) -> None:
    """bwa index"""
    runTool("bwa", ["bwa", "index", fasta, "-p", output_name])


def bwa(index: str, f1: str, f2: str, output_name: str, threads: int = 1) -> None:
    """Run bwa"""
    runTool(
        "bwa",
        [
            "bwa",
            "mem",
            "-t",
            str(threads),
            index,
            "-M",
            "-K",
            "100000000",
            "-v",
            "3",
            f1,
            f2,
            "-a",
            "-o",
            f"{output_name}.sam",
        ],
    )
    samtobam(output_name)


def extractDiploidCoverage(input_name: str, diploid_gene: str) -> str:
    region = regions_of_diploid[diploid_gene]
    output_name = input_name + f".diploid_gene_{diploid_gene}"
    logger.info(f"[WGS] Extract {diploid_gene} region to {output_name}")
    runTool(
        "samtools",
        [
            "samtools",
            "view",
            "-@",
            str(getThreads()),
            "-h",
            "-F",
            "1024",
            f"{input_name}.bam",
            region,
            "-o",
            f"{output_name}.bam",
        ],
    )

    # Calculate depth and save to TSV
    depth_name = output_name + ".depth"
    bam2Depth(f"{output_name}.bam", depth_name + ".tsv", get_all=False)

    # Calculate mean and std of depth
    df = readSamtoolsDepth(depth_name + ".tsv")
    mean = float(df["depth"].mean())
    std = float(df["depth"].std())

    # Save depth stat to JSON
    depth_stat_name = depth_name + ".stat"
    with open(depth_stat_name + ".json", "w") as f:
        json.dump(
            {"mean": mean, "std": std, "gene": diploid_gene, "name": input_name}, f
        )
    return depth_stat_name


def extractFromHg19(
    input_bam: str, output_name: str, hg19_type: str = "hs37d5", threads: int = 1
) -> None:
    """Extract records in KIR regions from hg19"""
    if hg19_type == "hs37d5":
        main_regions = "19:55200000-55400000 GL000209.1"
    elif hg19_type == "hs19":
        main_regions = "chr19:55200000-55400000 chr19_gl000209_random"
    else:
        raise NotImplementedError(
            "Our bam extraction only support hg19 and hs37d5 index"
        )

    logger.info(f"[WGS] Extract {main_regions} from {input_bam}")
    runTool(
        "samtools",
        ["samtools", "view", f"-@{threads}", "-h", "-F", "1024", input_bam]
        + main_regions.split()
        + ["-o", f"{output_name}.sam"],
    )
    samtobam(output_name)


def bam2fastq(
    input_bam: str, output_name: str, threads: int = 1, gzip: bool = True
) -> None:
    """
    Bam to fastq.

    Return format:
        * output_name.read.1.fq.gz
        * output_name.read.2.fq.gz
    """
    tmp_bam = output_name + ".sortn.bam"
    read_suffix = ".read.{}.fq"
    if gzip:
        read_suffix += ".gz"
    runTool(
        "samtools",
        ["samtools", "sort", "-@", str(threads), "-n", input_bam, "-o", tmp_bam],
    )
    runTool(
        "samtools",
        [
            "samtools",
            "fastq",
            "-@",
            str(threads),
            "-n",
            tmp_bam,
            "-1",
            f"{output_name}{read_suffix.format(1)}",
            "-2",
            f"{output_name}{read_suffix.format(2)}",
            "-0",
            "/dev/null",
            "-s",
            "/dev/null",
        ],
    )
    runShell(["rm", tmp_bam])
