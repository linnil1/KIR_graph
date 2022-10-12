import os
from glob import glob
from typing import Iterator, TextIO
from pathlib import Path
from functools import partial
from collections import defaultdict
from Bio.Seq import Seq
import pysam

from namepipe import nt, NameTask, compose, ConcurrentTaskExecutor
from kg_utils import runShell, threads, runDocker, samtobam
from graphkir.hisat2 import hisatMap


def downloadHG38(folder="index5"):
    output_name = f"{folder}/hg38.ucsc"
    if os.path.exists(f"{output_name}.fa.gz"):
        return output_name
    runShell("wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz "
             f"-O {output_name}.fa.gz")
    runShell("wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz "
             f"-O {output_name}.gff.gz")
    return output_name


def bwaIndex(input_name):
    output_name = input_name + ".bwa"
    if Path(output_name + ".bwt").exists():
        return output_name
    runDocker("bwa", f"bwa index {input_name}.fa.gz -p {output_name}")
    return output_name


def bwa(input_name, index):
    suffix = "." + index.replace("/", "_")
    output_name = input_name + suffix
    if Path(f"{output_name}.bam").exists():
        return output_name

    # main
    f1, f2 = input_name + ".read.1.fq.gz", input_name + ".read.2.fq.gz"
    if not os.path.exists(f1):
        f1, f2 = input_name + ".read.1.fq", input_name + ".read.2.fq"
    runDocker("bwa",
              f"bwa mem -t {threads} {index} -K 100000000 -p -v 3 "
              f" {f1} {f2} -o {output_name}.sam")
    samtobam(output_name)
    return output_name


# from kg_main import hisatMapWrap
def hisatMapWrap(input_name, index):
    # 1 to 1
    output_name = input_name + "." + index.replace("/", "_")
    if Path(f"{output_name}.bam").exists():
        return output_name
    f1, f2 = input_name + ".read.1.fq.gz", input_name + ".read.2.fq.gz"
    if not os.path.exists(f1):
        f1, f2 = input_name + ".read.1.fq", input_name + ".read.2.fq"
    hisatMap(index, f1, f2, output_name + ".bam", threads=threads)
    return output_name


def getTwbkSample() -> dict[str, str]:
    """ Find TWBK fastq """
    samples = glob(twbk_path + "*/WGS/FASTQ/*")
    samples_dict = defaultdict(list)
    for i in samples:
        samples_dict[i.split("/")[-1]].append(i)

    # filter
    for sample_id, paths in samples_dict.items():
        samples_dict[sample_id] = [i for i in paths if "10811-02" in i]
        assert samples_dict[sample_id]

    # from pprint import pprint
    # pprint(samples_dict)
    return {i: j[0] for i, j in samples_dict.items()}


def getFastQ(path: str) -> tuple[str, str] | None:
    """ Find fastq data """
    id = path.split('/')[-1]
    fastq = (
        f"{path}/{id}_S99_L999_R1_001.fastq.gz",
        f"{path}/{id}_S99_L999_R2_001.fastq.gz"
    )
    if not os.path.exists(fastq[0]) or not os.path.exists(fastq[1]):
        return None
    return fastq


def link(input_name: str, folder: str) -> str:
    output_name = folder + "/twbb"
    if os.path.exists(f"{output_name}.{input_name}.read.2.fq.gz"):
        return output_name + ".{}"

    global twbk_samples
    try:
        twbk_samples
    except NameError:
        twbk_samples = getTwbkSample()
    os.makedirs(folder, exist_ok=True)

    fastq = getFastQ(twbk_samples[input_name])
    if not fastq:
        print(f"Cannot found {input_name}")
        return output_name + ".{}"
    runShell(f"ln -s {fastq[0]} {output_name}.{input_name}.read.1.fq.gz")
    runShell(f"ln -s {fastq[1]} {output_name}.{input_name}.read.2.fq.gz")
    return output_name + ".{}"


def linkBam(input_name: str, folder: str) -> str:
    output_name = folder + "/hg19.twbb"
    if os.path.exists(f"{output_name}.{input_name}.bam.bai"):
        return output_name + ".{}"

    global twbk_bam_samples
    try:
        twbk_bam_samples
    except NameError:
        samples = glob(twbk_path + "*/WGS/hg19/BAM/GATK/*")
        samples_dict = defaultdict(list)
        for i in samples:
            name = i.split("/")[-1]
            samples_dict[name].append(f"{i}/{name}.hg19.sorted.realigned.maskDuplicates.recal")

        for name, paths in list(samples_dict.items()):
            samples_dict[name] = [i for i in paths if "10811-02" in i]

        twbk_bam_samples = samples_dict

    os.makedirs(folder, exist_ok=True)
    bam_path = twbk_bam_samples[input_name][0]
    if not os.path.exists(bam_path + ".bam"):
        print(input_name, "doesn't exists")
        return output_name + ".{}"

    runShell(f"ln -s {bam_path}.bam {output_name}.{input_name}.bam")
    runShell(f"ln -s {bam_path}.bam.bai {output_name}.{input_name}.bam.bai")
    return output_name + ".{}"


def linkHg38Bam(input_name: str, folder: str) -> str:
    output_name = folder + "/hg38.twbb"
    if os.path.exists(f"{output_name}.{input_name}.bam.bai"):
        return output_name + ".{}"

    # get input_name bam file
    global twbk_bam_samples
    try:
        twbk_bam_samples
    except NameError:
        samples = glob(twbk_path + "TWBR11002-01/WGS/GRCh38/BAM/GATKv2/*")
        samples_dict = defaultdict(list)
        for i in samples:
            name = i.split("/")[-1]
            samples_dict[name].append(f"{i}/{name}.hs38DH.dedup.postalt.sorted")
        twbk_bam_samples = samples_dict
    if input_name not in twbk_bam_samples:
        print(input_name, "doesn't exists")
        return output_name + ".{}"
    bam_path = twbk_bam_samples[input_name][0]

    # link
    os.makedirs(folder, exist_ok=True)
    runShell(f"ln -s {bam_path}.bam     {output_name}.{input_name}.bam")
    runShell(f"ln -s {bam_path}.bam.bai {output_name}.{input_name}.bam.bai")
    return output_name + ".{}"


def extractFromHg38(input_name):
    """ Extract hg38 """
    output_name = input_name + ".part.{}"
    output_all_name = input_name + ".part_merge"
    if os.path.exists(f"{output_all_name}.bam"):
        return output_all_name

    # samtools view -H /staging/biodata/lions/twbk/TWBR11002-01/WGS/GRCh38/BAM/GATKv2/NGS2_20150110G/NGS2_20150110G.hs38DH.dedup.postalt.sorted.bam | grep "SN:"| grep chr19 | awk '{print $2}' > hg38_chr19
    regions = ["chr19:50000000"] \
              + list(filter(None, map(lambda i: i.strip(), open("./hg38_chrun")))) \
              + list(filter(None, map(lambda i: i.strip(), open("./hg38_chr19"))))
    regions_text = " ".join(regions)
    b1 = f"{output_name.format('mapped')}.bam"
    b2 = f"{output_name.format('bad_pair')}.bam"
    runDocker("samtools", f"samtools view  -@ {threads} {input_name}.bam {regions_text} -o {b1}")
    # non proper pair but execlude depulicated
    runDocker("samtools", f"samtools view  -@ {threads} {input_name}.bam -F 2 -f 1024 -o {b2}")
    runDocker("samtools", f"samtools merge -@ {threads} -o {output_all_name}.sam {b1} {b2}")
    samtobam(output_all_name)
    return output_all_name


def extractFromHg19(input_name):
    """ Extract hg19 """
    output_name = input_name + ".part.{}"
    output_all_name = input_name + ".part_merge"
    if os.path.exists(f"{output_all_name}.bam"):
        return output_all_name

    regions = ("chr19:50000000 chr19_gl000208_random chr19_gl000209_random "
               "chrUn_gl000211 chrUn_gl000212 chrUn_gl000213 chrUn_gl000214 "
               "chrUn_gl000215 chrUn_gl000216 chrUn_gl000217 chrUn_gl000218 "
               "chrUn_gl000219 chrUn_gl000220 chrUn_gl000221 chrUn_gl000222 "
               "chrUn_gl000223 chrUn_gl000224 chrUn_gl000225 chrUn_gl000226 "
               "chrUn_gl000227 chrUn_gl000228 chrUn_gl000229 chrUn_gl000230 "
               "chrUn_gl000231 chrUn_gl000232 chrUn_gl000233 chrUn_gl000234 "
               "chrUn_gl000235 chrUn_gl000236 chrUn_gl000237 chrUn_gl000238 "
               "chrUn_gl000239 chrUn_gl000240 chrUn_gl000241 chrUn_gl000242 "
               "chrUn_gl000243 chrUn_gl000244 chrUn_gl000245 chrUn_gl000246 "
               "chrUn_gl000247 chrUn_gl000248 chrUn_gl000249 ")
    b1 = f"{output_name.format('mapped')}.bam"
    b2 = f"{output_name.format('bad_pair')}.bam"
    runDocker("samtools", f"samtools view  -@ {threads} {input_name}.bam {regions} -o {b1}")
    # non proper pair but execlude depulicated
    runDocker("samtools", f"samtools view  -@ {threads} {input_name}.bam -F 2 -f 1024 -o {b2}")
    runDocker("samtools", f"samtools merge -@ {threads} -o {output_all_name}.sam {b1} {b2}")
    samtobam(output_all_name)
    return output_all_name


def bam2fastqViaSamtools(input_name):
    """ Bam to fastq """
    output_name = input_name
    if os.path.exists(f"{output_name}.read.2.fq"):
        return output_name
    bam = input_name + ".sortn"
    runDocker("samtools", f"samtools sort -@ {threads} -n {input_name}.bam -o {bam}.bam")
    runDocker("samtools", f"samtools fastq -@ {threads} -n {bam}.bam "
                          f"-1 {output_name}.read.1.fq "
                          f"-2 {output_name}.read.2.fq "
                          f"-0 /dev/null -s /dev/null")
    return output_name


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
            if read_l.is_mapped:
                ref_l_name, ref_l_pos = read_l.reference_name, str(read_l.reference_start)
            else:
                ref_l_name, ref_l_pos = "*", "*"
            if read_r.is_mapped:
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
        return output_name
    bam2fastq(input_name + ".bam",
              output_name + ".read.1.fq",
              output_name + ".read.2.fq")
    return output_name


def trimHisatMap(input_name):
    """ Trim very large bam """
    output_name = input_name + ".trim"
    if os.path.exists(f"{output_name}.bam"):
        return output_name
    runDocker("samtools", f"samtools view -h -@ {threads} -F 4 {input_name}.bam -o {output_name}.sam")
    samtobam(output_name)
    return output_name


def copyBam(input_name, folder):
    """ Copy bamfile to persistence storage """
    os.makedirs(folder, exist_ok=True)
    runShell(f"cp {input_name}.bam {folder}/")
    runShell(f"cp {input_name}.bam.bai {folder}/")


if __name__ == "__main__":
    twbk_path = "/staging/biodata/lions/twbk/"
    ngs_sample = """
NGS2_20150412B
NGS2015038F
NGS2_20150405B
NGS2_20150211F
NGS2_20150110G
NGS2_20150505B
NGS20140612E
NGS20140704E
NGS2_20150503F
NGS2_20150603A
NGS2015036H
NGS2_20150207B
NGS20140605B
NGS20140609B
NGS20140602B
NGS2_20150110C
NGS2_20150406B
NGS2_20150406E
NGS2_20150407H
NGS2_20150105D
NGS2_20150107C
NGS2_20150304C
NGS2_20150603C
NGS2_20150112G
NGS2_20150502A
NGS2_20150306F
NGS2_20150110B
NGS2_20150602G
NGS2_20150102A
NGS2015021H
NGS2_20150312E
NGS2015041A
NGS20140611B
NGS2_20150301B
NGS2_20150101B
NGS2_20150309E
NGS20140710H
NGS2_20150405C
NGS2_20150507C
NGS2_20150502H
NGS20140709D
NGS2_20150101A
NGS2015012C
NGS20140710C
NGS20140605A
NGS20140603E
NGS2_20150201B
NGS20150111C
NGS2_20150301D
NGS2_20150304H
NGS20140702A
NGS2_20150108C
NGS2_20150201D
NGS20140608H
NGS2_20150510B
NGS2_20150311G
NGS2_20150311A
NGS2_20150506G
NGS20140609D
NGS20140708D
NGS2_20150103H
NGS2_20150302F
NGS2015035F
NGS2015015C
NGS2015023B
NGS2_20150106D
NGS2_20150102G
NGS20140801F
NGS2_20150306A
NGS2_20150111F
NGS2_20150209H
NGS2_20150508A
NGS2_20150208A
NGS2015032E
NGS20140701H
NGS2015026F
NGS20140801E
NGS20140603G
NGS2_20150312B
NGS2_20150305E
NGS2_20150211E
NGS2_20150104E
    """
    ngs_sample = list(filter(None, map(lambda i: i.strip(), ngs_sample.split("\n"))))
    print(ngs_sample)
    print(len(ngs_sample))
    ngs_sample = ngs_sample[0:7]
    ref = "hg38_ucsc"

    threads = 7
    if ref == "hg38_ucsc":  # fastq -> ucsc bam
        index = compose([
            "index",
            downloadHG38,
            bwaIndex,
        ])
        for id in ngs_sample:
            samples = link(id, "data_tmp")
        samples = compose([samples, partial(bwa, index=str(index))])

    NameTask.default_executor = ConcurrentTaskExecutor()
    NameTask.default_executor.threads = 7
    threads = 3

    if ref in ["hg38", "hg19", "hg38_ucsc"]:
        # bam -> fastq
        for id in ngs_sample:
            if ref == "hg38":
                samples = linkHg38Bam(id, "data_tmp")
            elif ref == "hg19":
                samples = linkBam(id, "data_tmp")
        print(samples)
        samples = compose([
            samples,
            extractFromHg38 if ref == "hg38" else extractFromHg19,
            # bam2fastqViaSamtools,
            bam2fastqWrap,
        ])
    elif ref == "fastq":
        # fastq
        for id in ngs_sample:
            samples = link(id, "data_tmp")
    print(samples)

    # fastq -> kir bam
    ref_index = "index/kir_2100_withexon_ab_2dl1s1.leftalign.mut01"
    index = ref_index + ".graph"
    compose([
        samples,
        partial(hisatMapWrap, index=str(index)),
        trimHisatMap,
        partial(copyBam, folder="data"),
    ])
