from pathlib import Path
from functools import partial
import pandas as pd

from namepipe import NameTask, compose, ConcurrentTaskExecutor, BaseTaskExecutor
from graphkir.utils import runShell

from kg_wgs import (
    downloadHg19,
    extractFromWGS,
    bam2fastqWrap,
    bam2fastqViaSamtools,
    extractHg19Depth,
    downloadHg38,
)
from kg_utils import back, SlurmTaskExecutor, addSuffix
from kg_mapping import bwa, bwaIndex, hisatMapWrap, trimBam


def linkHPRCSample(input_folder):
    df = pd.read_csv("hprc.csv")
    fastq_folder = "/staging/biology/zxc898977/rawData/WGS_Ncbi"
    df = df[df["sample_id"].notna()]
    print(df)
    # df = df[df["t1k"] == 1]
    # assert len(df) == 26
    for sample in df.itertuples():
        print(sample.id)
        print(sample.sample_id)
        assert Path(f"{fastq_folder}/{sample.sample_id}_1.fastq.gz").exists()
        if Path(f"{input_folder}/hprc.{sample.id}.read.2.fq.gz").exists():
            continue
        runShell(f"ln -s {fastq_folder}/{sample.sample_id}_1.fastq.gz {input_folder}/hprc.{sample.id}.read.1.fq.gz")
        runShell(f"ln -s {fastq_folder}/{sample.sample_id}_2.fastq.gz {input_folder}/hprc.{sample.id}.read.2.fq.gz")
    return input_folder + "/hprc.{}.read"


twbb_samples_str = """
NGS2_20150412B NGS2015038F    NGS2_20150405B NGS2_20150211F NGS2_20150110G
NGS2_20150505B NGS20140612E   NGS20140704E   NGS2_20150503F NGS2_20150603A
NGS2015036H    NGS2_20150207B NGS20140605B   NGS20140609B   NGS20140602B
NGS2_20150110C NGS2_20150406B NGS2_20150406E NGS2_20150407H NGS2_20150105D
NGS2_20150107C NGS2_20150304C NGS2_20150603C NGS2_20150112G NGS2_20150502A
NGS2_20150306F NGS2_20150110B NGS2_20150602G NGS2_20150102A NGS2015021H
NGS2_20150312E NGS2015041A    NGS20140611B   NGS2_20150301B NGS2_20150101B
NGS2_20150309E NGS20140710H   NGS2_20150405C NGS2_20150507C NGS2_20150502H
NGS20140709D   NGS2_20150101A NGS2015012C    NGS20140710C   NGS20140605A
NGS20140603E   NGS2_20150201B NGS20150111C   NGS2_20150301D NGS2_20150304H
NGS20140702A   NGS2_20150108C NGS2_20150201D NGS20140608H   NGS2_20150510B
NGS2_20150311G NGS2_20150311A NGS2_20150506G NGS20140609D   NGS20140708D
NGS2_20150103H NGS2_20150302F NGS2015035F    NGS2015015C    NGS2015023B
NGS2_20150106D NGS2_20150102G NGS20140801F   NGS2_20150306A NGS2_20150111F
NGS2_20150209H NGS2_20150508A NGS2_20150208A NGS2015032E    NGS20140701H
NGS2015026F    NGS20140801E   NGS20140603G   NGS2_20150312B NGS2_20150305E
NGS2_20150211E NGS2_20150104E """
twbb_samples = list(map(lambda i: i.strip(), twbb_samples_str.split()))


def linkTWBBBamSample(input_folder):
    samples = twbb_samples[0:14]
    twbk_path = "/staging/biodata/lions/twbk/TWBR10811-02/WGS/hg19/BAM/GATK/{name}/{name}.hg19.sorted.realigned.maskDuplicates.recal.bam"
    for name in samples:
        if Path(f"{input_folder}/hg19.twbb.{name}.bam").exists():
            continue
        assert Path(twbk_path.format(name=name)).exists()
        runShell(f"ln -s {twbk_path.format(name=name)}     {input_folder}/hg19.twbb.{name}.bam")
        runShell(f"ln -s {twbk_path.format(name=name)}.bai {input_folder}/hg19.twbb.{name}.bam.bai")
    return input_folder + "/hg19.twbb.{}"


def linkTWBBFastqSample(input_folder):
    twbk_path = "/staging/biodata/lions/twbk/TWBR10811-02/WGS/FASTQ/{name}/{name}_S99_L999_R{strand}_001.fastq.gz"
    samples = twbb_samples[0:7]
    for name in samples:
        if Path(f"{input_folder}/twbb.{name}.read.1.fq.gz").exists():
            continue
        assert Path(twbk_path.format(name=name, strand=1)).exists()
        runShell(f"ln -s {twbk_path.format(name=name, strand=1)} {input_folder}/twbb.{name}.read.1.fq.gz")
        runShell(f"ln -s {twbk_path.format(name=name, strand=2)} {input_folder}/twbb.{name}.read.2.fq.gz")
    return input_folder + "/twbb.{}.read"


if __name__ == "__main__":
    index_folder = "index"
    data_folder = "data_tmp"
    cohort = "twbb_fastq"
    direct_on_kir = False
    ref_index = "index5/kir_2100_withexon_ab_2dl1s1.leftalign.mut01"
    index = "index/kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph"
    search_other_region = True
    TAIWANIA = True

    # these step in run on Taiwania HPC
    if TAIWANIA:
        if cohort == "twbb_bam":
            sample = compose([
                data_folder,
                linkTWBBBamSample,
            ])
            direct_on_kir = False
            wgs_type = "hg19"
        elif cohort == "twbb_fastq":
            sample = compose([
                data_folder,
                linkTWBBFastqSample,
                back,
            ])
        elif cohort == "hprc_fastq":
            sample = compose([
                data_folder,
                linkHPRCSample,
                back,
            ])
        else:
            raise ValueError(f"{cohort} not found")

        if not direct_on_kir:
            genome_index = compose([
                index_folder,
                # downloadHg19,
                downloadHg38,
                bwaIndex,
            ])
            wgs_type = "hg37d5"
            wgs_type = "hg38"

        exe       = SlurmTaskExecutor(threads_per_sample=1,  template_file="research/taiwania.template")
        exe_large = SlurmTaskExecutor(threads_per_sample=14, template_file="research/taiwania.14.template")
        NameTask.set_default_executor(exe)
        if not direct_on_kir:
            if "fastq" in cohort:
                sample = compose([
                    sample,
                    NameTask(partial(bwa, index=str(genome_index)), executor=exe_large),
                ])
                # compose([sample, extractHg19Depth])
            sample = compose([
                sample,
                partial(extractFromWGS, wgs_type=wgs_type, loose=search_other_region),
                bam2fastqWrap if search_other_region else bam2fastqViaSamtools,
                back,
            ])
        sample = compose([
            sample,
            NameTask(partial(hisatMapWrap, index=str(index)), executor=exe_large),
            partial(trimBam, remove=True),
        ])
    else:
        # these step in run on our server
        NameTask.set_default_executor(ConcurrentTaskExecutor(10))
        from kg_main import extractVariant, bam2DepthWrap, filterDepthWrap, cnPredict, plotCNWrap, kirTyping, mergeKirResult
        from kg_utils import compareResult
        # samples = "data_real/twbb.{}"
        # samples = "data_real/hprc.{}"
        # samples = "data_real/hg19.twbb.{}.part_merge.annot_read"
        # samples = "data_real/hg19.twbb.{}.part_strict"
        # samples = "data_real/hprc.{}.index_hs37d5.bwa.part_merge.annot_read"
        # samples = "data_real/twbb.{}.index_hs37d5.bwa.part_merge.annot_read"
        # samples = "data_real/hprc.{}.index_hs37d5.bwa.part_strict"
        # samples = "data_real/twbb.{}.index_hs37d5.bwa.part_strict"
        # samples = "data_real/hprc26.{}.index_hs37d5.bwa.part_strict"
        samples = "data_real/hprc.{}.index_hs37d5.bwa.part_strict"
        samples += ".index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim"
        NameTask.set_default_executor(ConcurrentTaskExecutor(threads=20))
        variant = compose([
            samples,
            partial(extractVariant, ref_index=str(ref_index)),
        ])

        cn = compose([
            variant,
            partial(addSuffix, suffix=".no_multi"),
            bam2DepthWrap,
            partial(filterDepthWrap, ref_index=str(ref_index)),
            NameTask(cnPredict)  # .set_depended(-1),
        ])
        # cn >> NameTask(partial(plotCNWrap, per_sample=False, show_depth=False), depended_pos=[-1])

        # allele typing
        sample_possible_ans = "hprc_summary"  # from kg_from_kelvin.py
        typing = compose([
            variant,
            partial(kirTyping, cn_input_name=cn, allele_method="pv"),  # pv_exonfirst_1
            NameTask(mergeKirResult, depended_pos=[0]),
            partial(compareResult, sample_name=sample_possible_ans),
        ])
        runShell("stty echo opost")
