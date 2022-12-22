from typing import Iterator, TextIO
from pathlib import Path
from Bio.Seq import Seq
import pysam

from graphkir.utils import (
    runShell,
    samtobam,
    getThreads,
)
from kg_utils import (
    runDocker
)


def downloadHg19(index_folder):
    output_name = f"{index_folder}/hs37d5"
    if Path(output_name + ".fa.gz").exists():
        return output_name
    runShell("wget https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz "
             f"-O {output_name}.fa.gz")
    return output_name


def bam2fastqViaSamtools(input_name):
    """ Bam to fastq """
    output_name = input_name
    if Path(f"{output_name}.read.2.fq").exists():
        return output_name + ".read"
    bam = input_name + ".sortn"
    runDocker("samtools", f"samtools sort  -@ {getThreads()} -n {input_name}.bam -o {bam}.bam")
    runDocker("samtools", f"samtools fastq -@ {getThreads()} -n {bam}.bam "
                          f"-1 {output_name}.read.1.fq "
                          f"-2 {output_name}.read.2.fq "
                          f"-0 /dev/null -s /dev/null")
    return output_name + ".read"


def extractHg19Depth(input_name):
    output_name = input_name + ".part_strict.depth"
    if Path(f"{output_name}.tsv").exists():
        return output_name

    regions = ("19:55200000-55400000 GL000209.1").split()
    for i, reg in enumerate(regions):
        runDocker("samtools",
                  f"samtools depth -@ {getThreads()} -aa -r {reg} "
                  f"{input_name}.bam -o {output_name}.{i}.tsv")
    runShell("cat " + " ".join(f"{output_name}.{i}.tsv" for i in range(len(regions))) + f" > {output_name}.tsv")
    runShell("rm "  + " ".join(f"{output_name}.{i}.tsv" for i in range(len(regions))))
    return output_name


def extractFromHg19(input_name, hg19_type: str = "hs37d5", loose: bool = False):
    """ Extract hg19 """
    if loose:
        output_name = input_name + ".part.{}"
        output_all_name = input_name + ".part_merge"
    else:
        output_all_name = input_name + ".part_strict"

    if Path(f"{output_all_name}.bam").exists():
        return output_all_name

    if hg19_type == "hs37d5":
        if loose:
            main_regions = ("19:55000000                  GL000209.1 ")
        else:
            main_regions = ("19:55200000-55400000         GL000209.1 ")
        other_region = ("GL000191.1 GL000192.1 GL000193.1 GL000194.1 GL000195.1 "
                        "GL000196.1 GL000197.1 GL000198.1 GL000199.1 GL000200.1 "
                        "GL000201.1 GL000202.1 GL000203.1 GL000204.1 GL000205.1 "
                        "GL000206.1 GL000207.1 GL000208.1            GL000210.1 "
                        "GL000211.1 GL000212.1 GL000213.1 GL000214.1 GL000215.1 "
                        "GL000216.1 GL000217.1 GL000218.1 GL000219.1 GL000220.1 "
                        "GL000221.1 GL000222.1 GL000223.1 GL000224.1 GL000225.1 "
                        "GL000226.1 GL000227.1 GL000228.1 GL000229.1 GL000230.1 "
                        "GL000231.1 GL000232.1 GL000233.1 GL000234.1 GL000235.1 "
                        "GL000236.1 GL000237.1 GL000238.1 GL000239.1 GL000240.1 "
                        "GL000241.1 GL000242.1 GL000243.1 GL000244.1 GL000245.1 "
                        "GL000246.1 GL000247.1 GL000248.1 GL000249.1 "
                        "hs37d5 NC_007605 ")

    elif hg19_type == "hg19":
        if loose:
            main_regions = ("chr19:55000000         chr19_gl000208_random chr19_gl000209_random ")
        else:
            main_regions = ("chr19:55200000-55400000                      chr19_gl000209_random ")
        other_region = ("chrUn_gl000211 chrUn_gl000212 chrUn_gl000213 chrUn_gl000214 chrUn_gl000215 "
                        "chrUn_gl000216 chrUn_gl000217 chrUn_gl000218 chrUn_gl000219 chrUn_gl000220 "
                        "chrUn_gl000221 chrUn_gl000222 chrUn_gl000223 chrUn_gl000224 chrUn_gl000225 "
                        "chrUn_gl000226 chrUn_gl000227 chrUn_gl000228 chrUn_gl000229 chrUn_gl000230 "
                        "chrUn_gl000231 chrUn_gl000232 chrUn_gl000233 chrUn_gl000234 chrUn_gl000235 "
                        "chrUn_gl000236 chrUn_gl000237 chrUn_gl000238 chrUn_gl000239 chrUn_gl000240 "
                        "chrUn_gl000241 chrUn_gl000242 chrUn_gl000243 chrUn_gl000244 chrUn_gl000245 "
                        "chrUn_gl000246 chrUn_gl000247 chrUn_gl000248 chrUn_gl000249 ")
    else:
        raise ValueError("hg19_type doesn't correct")

    threads = getThreads()
    if loose:
        b1 = f"{output_name.format('mapped')}.bam"
        b2 = f"{output_name.format('chrun')  }.bam"
        b3 = f"{output_name.format('bad_pair')}.bam"
        runDocker("samtools", f"samtools view  -@ {threads} {input_name}.bam -F 1024 {main_regions} -o {b1}")
        runDocker("samtools", f"samtools view  -@ {threads} {input_name}.bam -F 1024 {other_region} -o {b2}")
        runDocker("samtools", f"samtools view  -@ {threads} {input_name}.bam -F 1024 -F 2 -o {b3}")
        runDocker("samtools", f"samtools merge -@ {threads} -o {output_all_name}.sam {b1} {b2} {b3}")
        samtobam(output_all_name)
    else:
        runDocker("samtools", f"samtools view  -@ {threads} {input_name}.bam -F 1024 {main_regions} -h -o {output_all_name}.sam")
        samtobam(output_all_name)
    return output_all_name


def readPair(bam_file: str) -> Iterator[tuple[pysam.AlignedSegment, pysam.AlignedSegment]]:
    """
    Extract paired read from bam

    Args:
      bam_file: The bam file path

    Yields:
      left_record(str), right_record(str)
    """
    num_reads = 0
    num_pairs = 0
    reads: dict[tuple[str, bool], pysam.AlignedSegment] = {}
    done_read_ids = set()

    for record in pysam.AlignmentFile(bam_file, "rb").fetch(until_eof=True):
        read_id = str(record.query_name)
        flag = record.flag
        num_reads += 1
        if read_id in done_read_ids:
            continue
        if record.is_read1 == record.is_read2:
            print("Strange case", record)
            continue
        if not record.query_sequence:  # bwa only preserve one record
            continue

        next_id = (read_id, not record.is_read1)
        if next_id in reads:
            next_record = reads[next_id]
            del reads[next_id]
            done_read_ids.add(read_id)
            num_pairs += 1
            yield record, next_record
        # if not find -> save in temp
        else:
            reads[(read_id, record.is_read1)] = record

    # summary
    print("Reads:", num_reads, "Pairs:", num_pairs)


def writeFastq(f: TextIO, id: str, record: pysam.AlignedSegment):
    # double check
    seq = record.query_sequence
    assert seq
    assert record.query_qualities
    qlt = pysam.array_to_qualitystring(record.query_qualities)
    assert qlt
    assert len(seq) == len(qlt)

    if record.is_reverse:
        seq = str(Seq(seq).reverse_complement())
        qlt = qlt[::-1]

    # write to file
    if record.is_read1:
        strand = "1"
    else:
        strand = "2"
    f.write(
        f"@{id}/{strand}\n"
        f"{seq}\n"
        "+\n"
        f"{qlt}\n"
    )


def bam2fastq(bam_file: str, fastq1_file: str, fastq2_file: str):
    """Convert bam to fastq, the mapping info are stored in readID"""
    with (
        open(fastq1_file, "w") as f_l,
        open(fastq2_file, "w") as f_r,
    ):
        for read_l, read_r in readPair(bam_file):
            if read_r.is_read1:
                read_l, read_r = read_r, read_l
            if not (read_l.is_read1 and read_r.is_read2):
                print("Strange case")
                print(read_l, read_r)
                continue
            if read_l.is_mapped:  # type: ignore
                ref_l_name, ref_l_pos = read_l.reference_name, str(read_l.reference_start)
            else:
                ref_l_name, ref_l_pos = "*", "*"
            if read_r.is_mapped:  # type: ignore
                ref_r_name, ref_r_pos = read_r.reference_name, str(read_r.reference_start)
            else:
                ref_r_name, ref_r_pos = "*", "*"
            new_name = (f"{read_l.query_name}"
                        f"|{ref_l_name}|{ref_l_pos}"
                        f"|{ref_r_name}|{ref_r_pos}")
            writeFastq(f_l, new_name, read_l)
            writeFastq(f_r, new_name, read_r)


def bam2fastqWrap(input_name):
    """ bam2fastq with one input """
    output_name = input_name + ".annot_read"
    if Path(output_name + ".read.1.fq").exists():
        return output_name + ".read"
    bam2fastq(input_name + ".bam",
              output_name + ".read.1.fq",
              output_name + ".read.2.fq")
    return output_name + ".read"
