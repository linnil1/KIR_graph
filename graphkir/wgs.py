"""
WGS index/mapping part of graphkir
"""
from .utils import runShell, runDocker, samtobam


def downloadHg19(index_folder: str) -> str:
    """Download hs37d5"""
    output_name = f"{index_folder}/hs37d5.fa.gz"
    runShell(
        "wget https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
        f"  -O {output_name}"
    )
    return output_name


def bwaIndex(fasta: str, output_name: str) -> None:
    """bwa index"""
    runDocker("bwa", f"bwa index {fasta} -p {output_name}")


def bwa(index: str, f1: str, f2: str, output_name: str, threads: int = 1) -> None:
    """Run bwa"""
    runDocker(
        "bwa",
        f"bwa mem -t {threads} {index} -M -K 100000000 -v 3 "
        f" {f1} {f2} -a -o {output_name}.sam",
    )
    samtobam(output_name)


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

    runDocker(
        "samtools",
        f"samtools view -@{threads} -h -F 1024"
        f"  {input_bam} {main_regions} -o {output_name}.sam",
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
    runDocker("samtools", f"samtools sort  -@ {threads} -n {input_bam} -o {tmp_bam}")
    runDocker(
        "samtools",
        f"samtools fastq -@ {threads} -n {tmp_bam}"
        f"  -1 {output_name}{read_suffix.format(1)}"
        f"  -2 {output_name}{read_suffix.format(2)}"
        f"  -0 /dev/null -s /dev/null",
    )
    runShell(f"rm {tmp_bam}")
