"""
WGS index/mapping part of graphkir
"""

import json
from pathlib import Path

from .utils import getThreads, runShell, samtobam, logger, downloadFile
from .external_tools import runTool
from .samtools_utils import bam2Depth, readSamtoolsDepth

# Diploid regions for different reference genomes
# hs37d5 (hg19) diploid regions
regions_of_diploid_hg19 = {
    "VDR": "12:48235320-48298777",
    "RYR1": "19:38924331-39078204",
    "EGFR": "7:55086710-55279321",
}

# hg38 (GRCh38) diploid regions
regions_of_diploid_hg38 = {
    "VDR": "chr12:47841537-47904994",
    "RYR1": "chr19:38433691-38587564",
    "EGFR": "chr7:55019017-55211628",
}

# Combined lookup by reference type
regions_of_diploid = {
    "hg19": regions_of_diploid_hg19,
    "hg38": regions_of_diploid_hg38,
}

def downloadHg19(index_folder: str) -> str:
    """Download hs37d5 (hg19)"""
    output_name = f"{index_folder}/hs37d5.fa.gz"
    logger.info(f"[WGS] Download {output_name}")
    url = "https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
    downloadFile(url, output_name)
    return output_name


def downloadHg38(index_folder: str) -> str:
    """Download GRCh38 no_alt analysis set (hg38)"""
    output_name = f"{index_folder}/hs38noalt.fa.gz"
    logger.info(f"[WGS] Download {output_name}")
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
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


def extractDiploidCoverage(input_name: str, diploid_gene: str, ref_type: str = "hg19") -> str:
    """Extract diploid gene region coverage for CN normalization.
    
    Args:
        input_name: Input BAM file prefix (without .bam extension)
        diploid_gene: Diploid gene name (VDR, RYR1, or EGFR)
        ref_type: Reference genome type (hg19 or hg38)
    
    Returns:
        Path to depth stat file prefix
    """
    # Normalize ref_type to lookup key
    if ref_type != "hg19" and ref_type != "hg38":
        raise ValueError(f"Unsupported reference type: {ref_type}")
    
    region = regions_of_diploid[ref_type][diploid_gene]
    output_name = input_name + f".diploid_gene_{diploid_gene}"
    logger.info(f"[WGS] Extract {diploid_gene} region ({region}) to {output_name}")
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


def extractKirRegion(
    input_bam: str, output_name: str, ref_type: str = "hg19", threads: int = 1
) -> None:
    """Extract records in KIR regions from reference genome."""
    
    # Define KIR regions for different reference genomes
    if ref_type == "hg19":
        # hg19/GRCh37: KIR region on chr19 + unplaced contig
        main_regions = "19:55200000-55400000 GL000209.1"
    elif ref_type == "hg38":
        # hg38/GRCh38: KIR region on chr19
        main_regions = "chr19:54720000-54870000"
    else:
        raise NotImplementedError(
            f"Reference type '{ref_type}' not supported. "
            "Supported types: hg19, hg38"
        )

    logger.info(f"[WGS] Extract KIR region ({main_regions}) from {input_bam}")
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
