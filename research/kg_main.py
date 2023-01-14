import os
from glob import glob
from pathlib import Path
import pandas as pd
from functools import partial
from Bio import SeqIO

from namepipe import NameTask, compose, ConcurrentTaskExecutor
from graphkir import kir_msa
from graphkir.msa2hisat import buildHisatIndex, msa2HisatReference
from graphkir.hisat2 import extractVariantFromBam, readExons
from graphkir.kir_cn import bam2Depth, filterDepth, loadCN, predictSamplesCN
from graphkir.kir_typing import selectKirTypingModel
from graphkir.plot import showPlot, plotCN, plotGeneDepths
from graphkir.msa_leftalign import genemsaLeftAlign
from graphkir.utils import (
    runShell,
    samtobam,
    getThreads,
    mergeAllele,
    readFromMSAs,
    setThreads,
)

from vg import (
    msa2vcf, buildVG, mapVG,
    vgindex2VGTube, writeVGTubeSettings, startVGTube, setupVGTube
)
from kg_utils import (
    runDocker, compareResult, linkSamples, getAnswerFile, addSuffix, back
)
from kg_mapping import bowtie2, bowtie2Index, bwa, bwaIndex, hisatMapWrap
from kg_wgs import downloadHg19
from kg_create_data import createSamplesAllele, createSamplesReads
from kg_create_fake_intron import createFakeIntronSample
from kg_create_novel import addNovelFromMsaWrap, updateNovelAnswer
from kg_typing_novel import typingNovel
from kg_create_exonseq_only import (
    extractPairReadsOnceInBed,
    calcExonToBed,
    bam2fastq,
)


def createSamplesWithAnswer(input_name, N=10, seed=2022):
    # 0 -> 1
    # "linnil1_syn_30x_seed87/linnil1_syn_30x_seed87"
    name = input_name + "/" + input_name + f"_s{seed}"
    output_name = name + ".{}"
    if Path(f"{name}_summary.csv").exists():
        return output_name
    createSamplesAllele(N, basename=name, seed=seed)
    return output_name


def createFakeSamplesWithAnswer(input_name, N=10, seed=2022, fake_num=1):
    # 0 -> 1
    # "linnil1_syn_30x_seed87/linnil1_syn_30x_seed87"
    name = input_name + "/" + input_name + f"_fakeintron{fake_num}_s{seed}"
    output_name = name + ".{}"
    if Path(f"{name}_summary.csv").exists():
        return output_name
    createFakeIntronSample(N, basename=name, seed=seed)
    return output_name


def createSampleFastq(input_name, depth=30, seed=1031):
    # 1 -> 1
    output_base = input_name + f".{depth}x_s{seed}"
    output_name = output_base + ".read"
    if Path(f"{output_base}.sam").exists():
        return output_name
    seed += int(input_name.template_args[0])
    createSamplesReads(input_name + ".fa", output_base, depth=depth, seed=seed)
    return output_name


def buildKirMsaWrap(input_name, msa_type="ab_2dl1s1"):
    output_name = f"{input_name}/kir_2100_{msa_type}"
    if len(glob(output_name + ".KIR*")):
        return output_name
    kir_msa.buildKirMsa(msa_type, output_name, mergeMSA=mergeMSA, threads=getThreads())
    return output_name


def buildMsaWithCds(input_name, msa_type="ab_2dl1s1"):
    exon_name = f"{input_name}/kir_2100_withexon"
    output_name = exon_name + "_" + msa_type

    if not len(glob(exon_name + ".KIR*")):
        kir_msa.buildKirMsa("ab", exon_name, full_length_only=False)

    if len(glob(output_name + ".KIR*")):
        return output_name
    kir_msa.buildKirMsa(msa_type, output_name, input_msa_prefix=exon_name, mergeMSA=mergeMSA, threads=getThreads())
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

    buildHisatIndex(input_name, output_name, threads=getThreads())
    return output_name


def msa2HisatReferenceWrap(input_name):
    output_name = input_name + ".mut01"
    if Path(output_name + ".haplotype").exists():
        return output_name
    msa2HisatReference(input_name, output_name)
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
    per_gene = False
    cn_cluster = "CNgroup"  # kde CNgroup
    cn_select = "p75"
    assume_3DL3_diploid = True
    suffix_cn = f".{cn_select}.{cn_cluster}"
    # suffix_cn += "_smallg"  # if not -> set CNgroup's dev_decay = 1 and not b2
    suffix_cn += "_b2"
    cluster_dev = "0.08"
    cluster_method_kwargs = {}
    if cn_cluster == "CNgroup" and cluster_dev != "0.08":
        suffix_cn += "_dev" + cluster_dev
        cluster_method_kwargs['base_dev'] = float(cluster_dev)
    # b2
    cluster_method_kwargs['start_base'] = 2

    if per_gene:
        assert ".{}." in input_name
        suffix_cn += "_per_gene"
    if ".{}." not in input_name and assume_3DL3_diploid:
        suffix_cn += "_assume3DL3"
    if ".{}." in input_name:
        suffix_cn += ".cohort"
    output_name = input_name + suffix_cn

    if ".{}." not in input_name:
        # SELECT * FROM depths GROUP BY SAMPLE
        if Path(output_name + ".json").exists():
            return output_name
        try:
            predictSamplesCN([input_name  + ".tsv"],
                             [output_name + ".tsv"],
                             cluster_method=cn_cluster,
                             cluster_method_kwargs=cluster_method_kwargs,
                             select_mode=cn_select,
                             assume_3DL3_diploid=assume_3DL3_diploid,
                             save_cn_model_path=output_name + ".json")
        except AssertionError as e:
            # SKIP this sample becuase 3DL3 assumption is fail
            print("fail", input_name, e)
            return None
    else:  # cohort
        # SELECT * FROM depths
        output_name1 = input_name.replace_wildcard("_merge_cn") + suffix_cn
        if Path(output_name1 + ".json").exists():
            return output_name
        predictSamplesCN([name             + ".tsv" for name in input_name.get_input_names()],
                         [name + suffix_cn + ".tsv" for name in input_name.get_input_names()],
                         per_gene=per_gene,
                         cluster_method=cn_cluster,
                         cluster_method_kwargs=cluster_method_kwargs,
                         select_mode=cn_select,
                         assume_3DL3_diploid=False,
                         save_cn_model_path=output_name1 + ".json")
    return output_name


def kirTyping(input_name, cn_input_name, allele_method="pv"):
    # setup
    top_n = 600
    error_correction = True
    assert len(input_name.template_args) == 1
    id = input_name.template_args[0]
    cn_name = cn_input_name.output_name.template.format(id)
    output_name_template = cn_input_name.output_name + "." + allele_method
    output_name_template += ".compare_sum"
    if error_correction:
        output_name_template += ".var_errcorr"
    if top_n != 300:
        output_name_template += f".top{top_n}"
    output_name = output_name_template.format(id)

    if Path(output_name + ".tsv").exists():
        return output_name_template
    # debug
    # if "02" != input_name.template_args[0]:
    #     return output_name_template
    t = selectKirTypingModel(allele_method, input_name + ".json", top_n=top_n, variant_correction=error_correction)
    if not Path(cn_name + ".tsv").exists():
        return None
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


def mergeKirResult(input_name, add_id: bool = True):
    output_name = input_name.replace_wildcard("_merge")
    names = input_name.get_input_names()
    mergeAllele(
        [name + ".tsv" for name in names],
        output_name + ".tsv"
    )
    if add_id:
        df = pd.read_csv(output_name + ".tsv", sep="\t")
        df["id"] = [name.template_args[-1] for name in names]
        df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def plotCNWrap(input_name, per_sample=False, show_depth=True):
    figs = []

    output_name = input_name.replace_wildcard("_merge_cn")
    if not per_sample and Path(output_name + ".json").exists():
        figs.extend(plotCN(output_name + ".json"))
    else:
        for name in input_name.get_input_names():
            # hard-coded
            if show_depth:
                figs.extend(plotGeneDepths(name[:name.find(".depth") + 6] + ".tsv", title=name))
            figs.extend(plotCN(name + ".json"))
            # if len(figs) > 10:
            #     break
    showPlot(figs)


def realignBlock(file, method, threads: int = 1):
    """ rewrite realignBlock in graphkir/kir_msa.py """
    if Path(f"{file}.{method}.fa").exists():
        return file + f".{method}"
    if method == "clustalo":
        return kir_msa.clustalo(file, threads)
    elif method == "muscle":
        return kir_msa.muscle(file, threads)
    else:
        raise NotImplementedError


def mergeMSA(genes: kir_msa.GenesMsa,
             method: str = "clustalo",
             tmp_prefix: str = "tmp",
             threads: int = 1) -> kir_msa.GenesMsa:
    """ Rewrite in graphkir/kir_msa.py """
    blocks = kir_msa.splitMsaToBlocks(genes)
    files = kir_msa.blockToFile(blocks, tmp_prefix=tmp_prefix)
    task = tmp_prefix + ".{}" >> NameTask(partial(realignBlock, method=method, threads=threads))
    files_name = list(task.output_name.get_input_names())
    files = {i.template_args[-1]: i for i in files_name}
    # original method
    # files = kir_msa.realignBlock(files, method)
    print(files)
    blocks = kir_msa.fileToBlock(files)
    msa = kir_msa.mergeBlockToMsa(blocks)
    kir_msa.isEqualMsa(genes, msa)
    return msa


def samtobamWrap(input_name, keep=False):
    if Path(f"{input_name}.bam").exists():
        return input_name
    samtobam(input_name, keep=keep)
    return input_name


def getMSABackbone(input_name):
    output_name = input_name + ".backbone"
    if Path(output_name + ".fa").exists():
        return output_name

    genes = readFromMSAs(str(input_name))
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

    genes = readFromMSAs(str(input_name))
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


def linkAnswer(input_name, old_name):
    # input: all to 1
    # return 1 to 1
    output_path = input_name.replace_wildcard("_summary")
    answer_path = old_name.replace_wildcard("_summary")
    if Path(output_path + ".csv").exists():
        return input_name
    relative = "../" * (len(Path(answer_path).parents) - 1)
    runShell(f"ln -s {relative}{answer_path}.csv {output_path}.csv")
    return input_name


def typingNovelWrap(input_name, msa_name, variant_name):
    output_name = input_name + ".novel"
    if "00" not in input_name.template_args:
        return output_name
    # if Path(output_name + ".fa").exists():
    #     return output_name
    typingNovel(
        variant_name=variant_name.replace("{}", input_name.template_args[-1]),
        msa_name=msa_name.replace("{}", input_name.template_args[-1]),
        result_name=input_name,
        output_name=output_name,
    )
    return output_name


if __name__ == "__main__":
    # setThreads(20)
    NameTask.set_default_executor(ConcurrentTaskExecutor(threads=10))
    index_folder = "index"
    extract_exon = False
    add_novel = False
    answer_folder = "linnil1_syn"
    data_folder = "data"

    cohort = "100"
    if cohort == "10":
        N = 10
        seed1 = 44
        seed2 = 444
        depth = 30
    elif cohort == "100":
        N = 100
        seed1 = 2022
        seed2 = 1031
        depth = 30
    elif cohort == "10fake":
        data_folder = "data5"
        N = 10
        seed1 = 1214
        seed2 = 2022
        depth = 30
    elif cohort == "10_50x":
        N = 10
        seed1 = 44
        seed2 = 444
        depth = 50

    msa_index = compose([
        index_folder,
        # partial(buildKirMsaWrap, msa_type="ab_2dl1s1"),  # merge split ab ab_2dl1s1
        partial(buildMsaWithCds, msa_type="ab_2dl1s1"),
        leftAlignWrap,
    ])

    Path(answer_folder).mkdir(exist_ok=True)
    if cohort == "10fake":
        samples = compose([
            answer_folder,
            partial(createFakeSamplesWithAnswer, N=N, seed=seed1),
        ])
    else:
        samples = compose([
            answer_folder,
            partial(createSamplesWithAnswer, N=N, seed=seed1),
        ])

    # novel
    if add_novel:
        # bed = samples >> calcExonToBed
        samples = compose([
            samples,
            partial(addNovelFromMsaWrap, msa_name=str(msa_index)),
            # partial(addNovelWrap, bed_name=bed.output_name),
            NameTask(partial(updateNovelAnswer, old_name=samples.output_name), depended_pos=[-1]),
        ])

    samples = compose([
        samples,
        partial(createSampleFastq, depth=depth, seed=seed2),
        back,
        NameTask(partial(linkAnswer, old_name=samples.output_name), depended_pos=[-1]),
    ])

    if extract_exon:
        bed = samples >> calcExonToBed
        samples = compose([
            samples,
            partial(samtobamWrap, keep=True),
            partial(extractPairReadsOnceInBed, bed_name=bed.output_name),
            bam2fastq,
            back,
            partial(linkSamples, data_folder=answer_folder, new_name=Path(samples.output_name).name + "_exon"),
            NameTask(partial(linkAnswer, old_name=samples.output_name), depended_pos=[-1]),
        ])

    samples_ori = samples
    samples = samples >> partial(linkSamples, data_folder=data_folder)

    """
    # check where the simulated reads map on hg19
    genome_index = compose([
        index_folder,
        downloadHg19,
        bwaIndex,
    ])
    compose([
        samples,
        partial(bwa, index=str(genome_index)),
    ])
    exit()
    """

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
        NameTask(cnPredict)  # .set_depended(-1),  # parameters written in cnPredict funcion
    ])
    cn >> NameTask(partial(plotCNWrap, per_sample=True, show_depth=False)).set_depended(0)
    exit()

    typing = compose([
        variant,
        partial(kirTyping, cn_input_name=cn, allele_method="pv"),  # pv pv_exonfirst_1 pv_exonfirst_1.2
    ])

    novel_funcs = [typing]
    if False:
        novel_funcs.append(
            partial(typingNovelWrap,
                    msa_name=msa_index.output_name,
                    variant_name=variant.output_name)
        )
    novel = compose(novel_funcs)

    compose([
        novel,
        NameTask(mergeKirResult, depended_pos=[0]),
        partial(compareResult, sample_name=samples_ori.output_name, input_fasta_name=novel.output_name, plot=False),
    ])

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
    bowtie2_index = backbone_index >> bowtie2Index
    mapping = compose([
        samples,
        partial(bowtie2, index=str(bowtie2_index)),
    ])
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
