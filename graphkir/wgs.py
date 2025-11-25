"""
WGS index/mapping part of graphkir
"""
from .utils import runShell, samtobam, logger, downloadFile
from .external_tools import runTool


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


def extractDiploidCoverage(input_name: str, diploid_gene: str) -> None:
    print(f"\nThe diploid gene is {diploid_gene}.\n")
    regions = {
        "VDR": "12:48235320-48298777",
        "RYR1": "19:38924331-39078204",
        "EGFR": "7:55086710-55279321",
    }
    region = regions[diploid_gene]
    print(f"region is {region}.\n")
    runTool(
        "samtools",
        [
            "samtools",
            "view",
            "-@8",
            "-h",
            "-F",
            "1024",
            f"{input_name}.bam",
            region,
            "-o",
            f"{input_name}.diploid_gene.sam",
        ],
    )
    samtobam(f"{input_name}.diploid_gene")
    output_name = input_name + ".diploid_info.txt"
    runTool(
        "samtools",
        [
            "sh",
            "-c",
            f"samtools depth {input_name}.diploid_gene.bam | awk '{{sum+=$3; sumsq+=$3*$3}} END {{print sum/NR; print sqrt(sumsq/NR - (sum/NR)**2)}}' > {output_name}",
        ],
    )
    runShell(
        ["rm", f"{input_name}.diploid_gene.bam", f"{input_name}.diploid_gene.bam.bai"]
    )


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
