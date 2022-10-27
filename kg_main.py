import os
from glob import glob
from pathlib import Path
import pandas as pd
from namepipe import NameTask, compose, ConcurrentTaskExecutor
from functools import partial
from Bio import SeqIO

from graphkir import kir_msa
from graphkir.msa2hisat import buildHisatIndex, msa2HisatReference
from graphkir.hisat2 import hisatMap, extractVariantFromBam, readExons
from graphkir.kir_cn import bam2Depth, filterDepth, loadCN, predictSamplesCN
from graphkir.kir_typing import selectKirTypingModel
from graphkir.main import mergeAllele
from graphkir.plot import showPlot, plotCN, plotGeneDepths
from graphkir.msa_leftalign import genemsaLeftAlign
from graphkir.utils import (
    runShell,
    samtobam,
    threads,
)

from vg import (
    msa2vcf, buildVG, mapVG,
    vgindex2VGTube, writeVGTubeSettings, startVGTube, setupVGTube
)
from kg_utils import runDocker
from kg_create_data import createSamplesAllele, createSamplesReads
from kg_eval import compareCohort, readPredictResult, readAnswerAllele
from kg_extract_exon_seq import (
    extractPairReadsOnceInBed,
    calcExonToBed,
    bam2fastq,
)


def createSamplesWithAnswer(input_name, N=10):
    # 0 -> 1
    # "linnil1_syn_30x_seed87/linnil1_syn_30x_seed87"
    name = input_name + "/" + input_name
    output_name = name + ".{}"
    if Path(f"{name}_summary.csv").exists():
        return output_name
    createSamplesAllele(N, basename=name, seed=878)
    return output_name


def createSampleFastq(input_name, depth=30):
    # 1 -> 1
    output_name = input_name + ".read"
    if Path(f"{input_name}.sam").exists():
        return output_name
    createSamplesReads(input_name, depth=depth, seed=878)
    return output_name


def linkSamples(input_name, data_folder, new_name=None, fastq=True, fasta=False, sam=False):
    # 1 -> 1
    # only work for unix relative folder
    if new_name is None:
        name = Path(input_name.template).name
    else:
        name = new_name
    output_name = str(Path(data_folder) / name)
    new_name = output_name.format(input_name.template_args[0])
    relative = "../" * (len(Path(new_name).parents) - 1)
    if fastq:
        if Path(new_name + ".read.1.fq").exists():
            return output_name
        runShell(f"ln -s {relative}{input_name}.read.1.fq {new_name}.read.1.fq")  # add .read is for previous naming
        runShell(f"ln -s {relative}{input_name}.read.2.fq {new_name}.read.2.fq")
    if fasta:
        if Path(new_name + ".fa").exists():
            return output_name
        runShell(f"ln -s {relative}{input_name}.fa {new_name}.fa")
    if sam:
        if Path(new_name + ".sam").exists():
            return output_name
        if Path(f"{input_name}.read..sam").exists():
            runShell(f"ln -s {relative}{input_name}.read..sam {new_name}.sam")
        else:
            runShell(f"ln -s {relative}{input_name}.sam {new_name}.sam")
    return output_name


def bowtie2Index(input_name):
    output_name = input_name + ".bowtie2"
    if Path(output_name + ".1.bt2").exists():
        return output_name
    runDocker("bowtie",
              f"bowtie2-build {input_name}.fa {output_name} --threads {threads}")
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
              f"bowtie2 {args} --threads {threads} -x {index} "
              f"-1 {f1} -2 {f2} -a -S {output_name}.sam")
    samtobam(output_name)
    return output_name


def bwaIndex(input_name="index/kir_2100_raw.mut01"):
    # brute force (remove mut01)
    output_name = input_name + ".bwa"
    if Path(output_name + ".bwt").exists():
        return output_name
    runDocker("bwa",
              f"bwa index {input_name}.fa -p {output_name}")
    return output_name


def bwa(input_name, index, use_arg="default"):
    suffix = "." + index.replace("/", "_")
    args = ""
    output_name = input_name + suffix
    if Path(f"{output_name}.bam").exists():
        return output_name

    # main
    f1, f2 = input_name + ".read.1.fq", input_name + ".read.2.fq"
    runDocker("bwa",
              f"bwa mem -t {threads} {index} "
              f" {f1} {f2} -a -o {output_name}.sam")
    samtobam(output_name)
    return output_name


def buildKirMsaWrap(input_name, msa_type="ab_2dl1s1"):
    output_name = f"{input_name}/kir_2100_{msa_type}"
    if len(glob(output_name + ".KIR*")):
        return output_name
    kir_msa.buildKirMsa(msa_type, output_name, mergeMSA=mergeMSA)
    return output_name


def buildMsaWithCds(input_name, msa_type="ab_2dl1s1"):
    exon_name = f"{input_name}/kir_2100_withexon"
    output_name = exon_name + "_" + msa_type

    if not len(glob(exon_name + ".KIR*")):
        kir_msa.buildKirMsa("ab", exon_name, full_length_only=False)

    if len(glob(output_name + ".KIR*")):
        return output_name
    kir_msa.buildKirMsa(msa_type, output_name, input_msa_prefix=exon_name, mergeMSA=mergeMSA)
    return output_name


def leftAlignWrap(input_name):
    output_name = input_name + ".leftalign"
    if len(glob(output_name + ".KIR*")):
        return output_name
    genemsaLeftAlign(input_name, output_name)
    return output_name


def buildHisatIndexWrap(input_name):
    output_name = input_name + ".graph"
    if Path(output_name + ".8.ht2").exists():
        return output_name

    buildHisatIndex(input_name, output_name)
    return output_name


def msa2HisatReferenceWrap(input_name):
    output_name = input_name + ".mut01"
    if Path(output_name + ".haplotype").exists():
        return output_name
    msa2HisatReference(input_name, output_name)
    return output_name


def hisatMapWrap(input_name, index):
    # 1 to 1
    output_name = input_name + "." + index.replace("/", "_")
    if Path(f"{output_name}.bam").exists():
        return output_name
    f1, f2 = input_name + ".read.1.fq.gz", input_name + ".read.2.fq.gz"
    if not Path(f1).exists():
        f1, f2 = input_name + ".read.1.fq", input_name + ".read.2.fq"
    hisatMap(index, f1, f2, output_name + ".bam", threads=threads)
    return output_name


def extractVariant(input_name, ref_index):
    error_correction = False
    output_name = input_name + ".variant"
    if not error_correction:
        output_name += ".noerrcorr"
    if Path(output_name + ".json").exists():
        return output_name
    extractVariantFromBam(ref_index, input_name + ".bam", output_name, error_correction=error_correction)
    return output_name


def bam2DepthWrap(input_name):
    output_name = input_name + ".depth"
    if Path(output_name + ".tsv").exists():
        return output_name
    bam2Depth(input_name + ".bam", output_name + ".tsv")
    return output_name


def filterDepthWrap(input_name, ref_index, exon=False):
    output_name = input_name
    exon_regions = {}
    if exon:
        exon_regions = readExons(ref_index)
        output_name += ".exon"
    if Path(output_name + ".tsv").exists():
        return output_name
    filterDepth(input_name + ".tsv",
                output_name + ".tsv",
                exon_regions)
    return output_name


def cnPredict(input_name):
    cn_cluster = "CNgroup"
    cn_select = "p75"
    assume_3DL3_diploid = True
    suffix_cn = f".{cn_select}.{cn_cluster}"
    if ".{}." not in input_name and assume_3DL3_diploid:
        suffix_cn += "_assume3DL3"
    if ".{}." in input_name:
        suffix_cn += ".cohort"
    output_name = input_name + suffix_cn

    if ".{}." not in input_name:
        if Path(output_name + ".json").exists():
            return output_name
        predictSamplesCN([input_name  + ".tsv"],
                         [output_name + ".tsv"],
                         cluster_method=cn_cluster,
                         select_mode=cn_select,
                         assume_3DL3_diploid=assume_3DL3_diploid,
                         save_cn_model_path=output_name + ".json")
    else:  # cohort
        output_name1 = input_name.replace_wildcard("_merge_cn") + suffix_cn
        if Path(output_name1 + ".json").exists():
            return output_name
        predictSamplesCN([name             + ".tsv" for name in input_name.get_input_names()],
                         [name + suffix_cn + ".tsv" for name in input_name.get_input_names()],
                         cluster_method=cn_cluster,
                         select_mode=cn_select,
                         assume_3DL3_diploid=False,
                         save_cn_model_path=output_name1 + ".json")
    return output_name


def kirTyping(input_name, cn_input_name, allele_method="pv"):
    # setup
    top_n = 300
    assert len(input_name.template_args) == 1
    id = input_name.template_args[0]
    cn_name = cn_input_name.output_name.template.format(id)
    output_name_template = cn_input_name.output_name + "." + allele_method
    output_name_template += ".compare_sum"
    if top_n != 300:
        output_name_template += f".top{top_n}"
    output_name = output_name_template.format(id)

    if Path(output_name + ".tsv").exists():
        return output_name_template
    # debug
    # if "02" != input_name.template_args[0]:
    #     return output_name_template

    t = selectKirTypingModel(allele_method, input_name + ".json", top_n=top_n)
    cn = loadCN(cn_name + ".tsv")
    called_alleles = t.typing(cn)
    t.save(output_name + ".json")

    print(input_name)
    print(called_alleles)
    pd.DataFrame([{
        'name': output_name,
        'alleles': "_".join(called_alleles),
    }]).to_csv(output_name + ".tsv", sep="\t", index=False)
    return output_name_template


def kirResult(input_name, answer):
    output_name = input_name.replace_wildcard("_merge")

    mergeAllele(
        [name + ".tsv" for name in input_name.get_input_names()],
        output_name + ".tsv"
    )
    if Path(f"{answer}/{answer}.summary.csv").exists():
        answer = readAnswerAllele(f"{answer}/{answer}.summary.csv")
    else:
        answer = readAnswerAllele(f"{answer}/{answer}_summary.csv")
    predit = readPredictResult(output_name + ".tsv")
    print(output_name + ".tsv")
    # compareCohort(answer, predit, skip_empty=True)
    compareCohort(answer, predit, skip_empty=True, plot=True)
    return output_name


def plotCNWrap(input_name):
    figs = []

    output_name = input_name.replace_wildcard("_merge_cn")
    if Path(output_name + ".json").exists():
        figs.extend(plotCN(output_name + ".json"))
    else:
        for name in input_name.get_input_names():
            # hard-coded
            figs.extend(plotGeneDepths(name[:name.find(".depth") + 6] + ".tsv", title=name))
            figs.extend(plotCN(name + ".json"))
    showPlot(figs)


def realignBlock(file, method):
    """ rewrite realignBlock in graphkir/kir_msa.py """
    if Path(f"{file}.{method}.fa").exists():
        return file + f".{method}"
    if method == "clustalo":
        return kir_msa.clustalo(file)
    elif method == "muscle":
        return kir_msa.muscle(file)
    else:
        raise NotImplementedError


def mergeMSA(genes: kir_msa.GenesMsa,
             method: str = "clustalo",
             tmp_prefix: str = "tmp") -> kir_msa.Genemsa:
    """ Rewrite in graphkir/kir_msa.py """
    blocks = kir_msa.splitMsaToBlocks(genes)
    files = kir_msa.blockToFile(blocks, tmp_prefix=tmp_prefix)
    task = tmp_prefix + ".{}" >> NameTask(partial(realignBlock, method=method))
    files_name = list(task.output_name.get_input_names())
    files = {i.template_args[-1]: i for i in files_name}
    # original method
    # files = kir_msa.realignBlock(files, method)
    print(files)
    blocks = kir_msa.fileToBlock(files)
    msa = kir_msa.mergeBlockToMsa(blocks)
    kir_msa.isEqualMsa(genes, msa)
    return msa


def samtobamWrap(input_name):
    if Path(f"{input_name}.bam").exists():
        return input_name
    return samtobam(input_name)


def back(input_name):
    return str(Path(input_name.template).with_suffix(""))


def addSuffix(input_name, suffix):
    return input_name + suffix


def getMSABackbone(input_name):
    output_name = input_name + ".backbone"
    if Path(output_name + ".fa").exists():
        return output_name

    genes = kir_msa.readFromMSAs(str(input_name))
    sequences = []
    for gene, msa in genes.items():
        ref_name = msa.get_reference()[0]
        assert "BACKBONE" in ref_name
        sequences.extend(
            msa.select_allele([ref_name]).to_records(gap=False)
        )
    SeqIO.write(sequences, output_name + ".fa", "fasta")
    return output_name


def getMSAFullSequence(input_name):
    output_name = input_name + ".full"
    if Path(output_name + ".fa").exists():
        return output_name

    genes = kir_msa.readFromMSAs(str(input_name))
    sequences = []
    for gene, msa in genes.items():
        ref_name = msa.get_reference()[0]
        assert "BACKBONE" in ref_name
        msa = msa.remove(ref_name)
        sequences.extend(
            msa.to_records(gap=False)
        )
    SeqIO.write(sequences, output_name + ".fa", "fasta")
    return output_name


if __name__ == "__main__":
    NameTask.set_default_executor(ConcurrentTaskExecutor(threads=20))

    data_folder = "data5"
    index_folder = "index5"
    Path(data_folder).mkdir(exist_ok=True)
    Path(index_folder).mkdir(exist_ok=True)
    extract_exon = False

    answer_folder = "linnil1_syn_20x"
    answer_folder = "linnil1_syn_30x_seed87"
    Path(answer_folder).mkdir(exist_ok=True)

    samples = compose([
        answer_folder,
        partial(createSamplesWithAnswer, N=10),
        partial(createSampleFastq, depth=20),
        back,
    ])

    if extract_exon:
        new_answer_folder = answer_folder + "_exon"
        Path(new_answer_folder).mkdir(exist_ok=True)
        samples = compose([
            samples,
            partial(linkSamples, data_folder=new_answer_folder, fastq=False, sam=True, fasta=True),
        ])
        bed = samples >> calcExonToBed
        samples = compose([
            samples,
            samtobamWrap,
            partial(extractPairReadsOnceInBed, bed_name=bed),
            bam2fastq,
            back,
            partial(linkSamples, data_folder=new_answer_folder, new_name=new_answer_folder + ".{}"),
        ])
        runShell(f"ln -fs ../{answer_folder}/{answer_folder}_summary.csv {new_answer_folder}/{new_answer_folder}_summary.csv")
        answer_folder = new_answer_folder

    samples = samples >> partial(linkSamples, data_folder=data_folder)

    msa_index = compose([
        index_folder,
        partial(buildKirMsaWrap, msa_type="ab_2dl1s1"),
        # partial(buildMsaWithCds, msa_type="ab_2dl1s1"),
        leftAlignWrap,
    ])
    ref_index = msa_index >> msa2HisatReferenceWrap
    index = ref_index >> buildHisatIndexWrap
    variant = compose([
        samples,
        partial(hisatMapWrap, index=str(index)),
        partial(extractVariant, ref_index=str(ref_index)),
    ])

    cn = compose([
        variant,
        partial(addSuffix, suffix=".no_multi"),
        bam2DepthWrap,
        partial(filterDepthWrap, ref_index=str(ref_index), exon=extract_exon),
        NameTask(cnPredict)  # .set_depended(0),
    ])
    typing = compose([
        variant,
        partial(kirTyping, cn_input_name=cn, allele_method="pv"),  # pv_exonfirst_1
    ])

    typing >> partial(typingNovelWrap,
        msa_name=msa_index.output_name,
        variant_name=variant.output_name,
    )

    typing >> NameTask(partial(kirResult, answer_name=samples_ori.output_name), depended_pos=[0])
    # cn >> nt(plotCNWrap).set_depended(0)

    """
    # vg
    vg_index = compose([
        msa_index,
        getMSABackbone,
        back,
        partial(addSuffix, suffix=".{}"),
        msa2vcf,
        NameTask(buildVG, depended_pos=[-1]),
    ])
    mapping = compose([
        samples,
        partial(mapVG, index=str(vg_index)),
    ])

    # visualize
    tube_index = vg_index >> vgindex2VGTube
    tube = compose([None, setupVGTube])
    samples_vis = compose([
        samples,
        partial(mapVG, index=str(vg_index), output_type="gam"),
        NameTask(partial(writeVGTubeSettings, index=str(tube_index), tube_path=str(tube)), depended_pos=[-1]),
        NameTask(partial(startVGTube, tube_path=str(tube), port=8001), depended_pos=[-1]),
    ])

    # bowtie or bwa
    backbone_index = compose([
        msa_index,
        getMSABackbone,
        # getMSAFullSequence,
    ])

    # bowtie
    mapping = compose([
        samples,
        partial(bowtie2, index=str(bowtie2_index)),
    ])
    bowtie2_index = backbone_index >> bowtie2Index
    # bowtie2_index = index >> bowtie2BuildConsensus  # "index/kir_2100_?_cons"
    # bowtie2_index = index >> "index/kir_2100_raw.mut01" >> bowtie2BuildFull >> "index/kir_2100_raw_full" # index = "index/kir_2100_raw_full"
    # mapping = samples >> bowtie2.set_args(index=str(bowtie2_index), use_arg="ping")

    # bwa
    bwa_index = backbone_index >> bwaIndex
    mapping = compose([
        samples,
        partial(bwa, index=str(bwa_index)),
    ])
    """
    runShell("stty echo opost")
