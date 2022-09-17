import os
from glob import glob
from pathlib import Path
import pandas as pd
from namepipe import nt, NameTask

from graphkir.gk_build_index import buildHisatIndex, msa2HisatReference
from graphkir.gk_build_msa import buildKirMsa
from graphkir.gk_hisat2 import hisatMap, extractVariantFromBam
from graphkir.gk_cn import predictSamplesCN, loadCN
from graphkir.gk_kir_typing import selectKirTypingModel
from graphkir.gk_main import mergeAllele
from graphkir.gk_plot import showPlot, plotCN
from graphkir.gk_msa_leftalign import genemsaLeftAlign
from graphkir.gk_utils import (
    runShell,
    samtobam,
)

from kg_utils import (
    threads,
    runDocker,
)
import kg_create_data
from kg_eval import compareCohort, readPredictResult, readAnswerAllele
from kg_extract_exon_seq import extractExonPairReads


@nt
def link10Samples(input_name):
    output_name = input_name.template.format("test10.{}")
    if not input_name.template_args[0].isdigit() or int(input_name.template_args[0]) >= 10:
        return output_name

    id = input_name.template_args[0]
    name = output_name.format(id)
    if Path(f"{name}.read.1.fq").exists():
        return output_name

    runShell(f"ln -s {Path(input_name).name}.read.1.fq {output_name.format(id)}.read.1.fq")
    runShell(f"ln -s {Path(input_name).name}.read.2.fq {output_name.format(id)}.read.2.fq")
    return output_name


@nt
def createSamples(input_name):
    # 0 -> 1
    # "linnil1_syn_30x_seed87/linnil1_syn_30x_seed87"
    name = input_name + "/" + input_name
    output_name = name + ".{}.read"
    N = 10
    for i in range(N):
        if Path(f"{name}.{i:02d}.read.1.fq").exists():
            return output_name
    kg_create_data.createSamples(N=10, basename=name, depth=30, seed=87)
    return output_name


@nt
def linkSamples(input_name, data_folder):
    # 1 -> 1
    name = input_name.split('/')[0]
    output_name = os.path.join(data_folder, name + ".{}")
    new_name = output_name.format(input_name.template_args[0])
    if Path(new_name + ".read.1.fq").exists():
        return output_name
    runShell(f"ln -s ../{input_name}.1.fq {new_name}.read.1.fq")  # add .read is for previous naming
    runShell(f"ln -s ../{input_name}.2.fq {new_name}.read.2.fq")
    return output_name


@nt
def bowtie2BuildFull(input_name):
    # No matter what index, this file will have the same sequences
    name = input_name
    assert "kir_2100_raw" in name
    if ".mut01" in name:
        name = name.replace(".mut01", "")
    name += "_full"
    if Path(name + ".1.bt2").exists():
        return name
    runDocker("bowtie",
              f"bowtie2-build {input_name}_sequences.fa {name} --threads {threads}")
    return name


@nt
def bowtie2BuildConsensus(input_name="index/kir_2100_raw.mut01"):
    # brute force (remove mut01)
    name = input_name
    if ".mut01" in name:
        name = name.replace(".mut01", "")
    name += "_cons"
    if Path(name + ".1.bt2").exists():
        return name
    runDocker("bowtie",
              f"bowtie2-build {input_name}_backbone.fa {name} --threads {threads}")
    return name


@nt
def bowtie2(input_name, index, use_arg="default"):
    suffix = "." + index.replace("/", "_") + ".bowtie2"
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


@nt
def bwaIndex(input_name="index/kir_2100_raw.mut01"):
    # brute force (remove mut01)
    name = input_name
    if ".mut01" in name:
        name = name.replace(".mut01", "")

    name += "_cons_bwa"
    if Path(name + ".bwt").exists():
        return name
    runDocker("bwa",
              f"bwa index {input_name}_backbone.fa -p {name}")
    return name


@nt
def bwa(input_name, index, use_arg="default"):
    suffix = "." + index.replace("/", "_") + ".bwa"
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


@nt
def buildKirMsaWrap(input_name, msa_type="ab_2dl1s1"):
    output_name = f"{input_name}/kir_2100_{msa_type}"
    if len(glob(output_name + ".KIR*")):
        return output_name
    buildKirMsa(msa_type, output_name)
    return output_name


@nt
def leftAlignWrap(input_name):
    output_name = input_name + ".leftalign"
    if len(glob(output_name + ".KIR*")):
        return output_name
    genemsaLeftAlign(input_name, output_name)
    return output_name


@nt
def buildHisatIndexWrap(input_name):
    output_name = input_name + ".graph"
    if Path(output_name + ".8.ht2").exists():
        return output_name

    buildHisatIndex(input_name, output_name)
    return output_name


@nt
def msa2HisatReferenceWrap(input_name):
    output_name = input_name + ".mut01"
    if Path(output_name + ".haplotype").exists():
        return output_name
    msa2HisatReference(input_name, output_name)
    return output_name


@nt
def hisatMapWrap(input_name, index):
    # 1 to 1
    output_name = input_name + "." + index.replace("/", "_")
    if Path(f"{output_name}.bam").exists():
        return output_name
    f1, f2 = input_name + ".read.1.fq", input_name + ".read.2.fq"
    hisatMap(index, f1, f2, output_name + ".bam")
    return output_name


@nt
def extractVariant(input_name, ref_index):
    output_name = input_name + ".variant"
    if Path(output_name + ".json").exists():
        return output_name
    extractVariantFromBam(ref_index, input_name + ".bam", output_name)
    return output_name


@nt
def cnPredict(input_name, ref_index, exon=False):
    suffix_variant = ".no_multi"
    cn_cluster = "CNgroup"
    cn_select = "p75"
    assume_3DL3_diploid = True

    suffix_depth = f"{suffix_variant}.depth{'.exon' if exon else ''}"
    suffix_cn = f"{suffix_depth}.{cn_select}.{cn_cluster}"
    if ".{}." not in input_name and assume_3DL3_diploid:
        suffix_cn += "_assume3DL3"
    if ".{}." in input_name:
        suffix_cn += ".cohort"
    output_name = input_name + suffix_cn
    exon_regions = {}
    if exon:
        exon_regions = readExons(ref_index)

    if ".{}." not in input_name:
        if Path(output_name + ".json").exists():
            return output_name
        predictSamplesCN([input_name + suffix_variant + ".bam"],
                         [input_name + suffix_cn      + ".tsv"],
                         bam_selected_regions=exon_regions,
                         cluster_method=cn_cluster,
                         select_mode=cn_select,
                         assume_3DL3_diploid=assume_3DL3_diploid,
                         save_cn_model_path=output_name + ".json")
    else:  # cohort
        output_name1 = input_name.replace_wildcard("_merge_cn") + suffix_cn
        if Path(output_name1 + ".json").exists():
            return output_name
        predictSamplesCN([name + suffix_variant + ".bam" for name in input_name.get_input_names()],
                         [name + suffix_cn      + ".tsv" for name in input_name.get_input_names()],
                         bam_selected_regions=exon_regions,
                         cluster_method=cn_cluster,
                         select_mode=cn_select,
                         assume_3DL3_diploid=False,
                         save_cn_model_path=output_name1 + ".json")
    return output_name


@nt
def kirTyping(input_name, cn_input_name, allele_method="pv"):
    # setup
    assert len(input_name.template_args) == 1
    id = input_name.template_args[0]
    cn_name = cn_input_name.output_name.template.format(id)
    output_name_template = cn_input_name.output_name + "." + allele_method
    output_name_template += ".compare_sum"
    output_name = output_name_template.format(id)

    if Path(output_name + ".tsv").exists():
        return output_name_template
    # debug
    # if "02" != input_name.template_args[0]:
    #     return output_name_template

    t = selectKirTypingModel(allele_method, input_name + ".json")
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


@nt
def kirResult(input_name, answer):
    output_name = input_name.replace_wildcard("_merge")

    mergeAllele(
        [name + ".tsv" for name in input_name.get_input_names()],
        output_name + ".tsv"
    )

    answer = readAnswerAllele(f"{answer}/{answer}.summary.csv")
    predit = readPredictResult(output_name + ".tsv")
    compareCohort(answer, predit, skip_empty=True)
    return output_name


@nt
def extractExon(input_name, folder):
    Path(folder).mkdir(exist_ok=True)
    copy_name = f"{folder}/{Path(input_name).name.replace('.read', '')}"
    output_template = f"{folder}/{Path(input_name.template).name.replace('.read', '').replace('.{}', '_exon.{}')}"
    output_name = output_template.format(input_name.template_args[0])

    if Path(f"{output_name}.read.1.fq").exists():
        return output_template + ".read"

    runShell(f"ln -s ../{input_name.replace('.read', '')}.fa {copy_name}.fa")
    runShell(f"ln -s ../{input_name}..sam                    {copy_name}.sam")
    extractExonPairReads(f"{copy_name}.fa", f"{copy_name}.sam", output_name)
    return output_template + ".read"


@nt
def plotCNWrap(input_name):
    figs = []

    output_name = input_name.replace_wildcard("_merge_cn")
    if Path(output_name + ".json").exists():
        figs.extend(plotCN(output_name + ".json"))
    else:
        for name in input_name.get_input_names():
            figs.extend(plotCN(name + ".json"))
    showPlot(figs)


if __name__ == "__main__":
    data_folder = "data5"
    index_folder = "index5"
    Path(data_folder).mkdir(exist_ok=True)
    Path(index_folder).mkdir(exist_ok=True)
    extract_exon = False

    answer_folder = "linnil1_syn_wide"
    answer_folder = "linnil1_syn_exon"
    answer_folder = "linnil1_syn_30x"
    answer_folder = "linnil1_syn_30x_seed87"
    Path(answer_folder).mkdir(exist_ok=True)

    samples = answer_folder >> createSamples

    if extract_exon:
        samples = samples >> extractExon.set_args(folder=answer_folder + "_exon")
        runShell(f"ln -fs ../{answer_folder}/{answer_folder}.summary.csv {answer_folder}_exon/{answer_folder}_exon.summary.csv")
        answer_folder = answer_folder + "_exon"

    samples = samples >> linkSamples.set_args(data_folder)
    if answer_folder == "linnil1_syn_wide":
        samples = samples >> link10Samples

    msa_index = index_folder >> buildKirMsaWrap.set_args("ab_2dl1s1") >> leftAlignWrap
    ref_index = msa_index >> msa2HisatReferenceWrap
    index = ref_index >> buildHisatIndexWrap

    print(index)
    # samples = "data2/linnil1_syn_30x_seed87.{}"
    mapping = samples >> hisatMapWrap.set_args(index=str(index))
    variant = mapping >> extractVariant.set_args(ref_index=str(ref_index))

    cn = variant >> cnPredict.set_args(ref_index=str(ref_index), exon=extract_exon)  # .set_depended(0)
    # cn = variant >> cnPredict.set_args(ref_index=str(ref_index), exon=extract_exon).set_depended(0)

    typing = variant >> kirTyping.set_args(cn, "pv_exonfirst") >> kirResult.set_args(answer=answer_folder).set_depended(0)

    # cn >> plotCNWrap.set_depended(0)

    # bowtie mapping rate
    # bowtie2_index = index >> bowtie2BuildConsensus  # "index/kir_2100_?_cons"
    # bowtie2_index = index >> "index/kir_2100_raw.mut01" >> bowtie2BuildFull >> "index/kir_2100_raw_full" # index = "index/kir_2100_raw_full"
    # mapping = samples >> bowtie2.set_args(index=str(bowtie2_index))
    # mapping = samples >> bowtie2.set_args(index=str(bowtie2_index), use_arg="ping")

    # bwa
    # bwa_index = index >> bwaIndex  # "index/kir_2100_?_cons"
    # mapping = samples >> bwa.set_args(index=str(bwa_index))
    runShell("stty echo opost")
