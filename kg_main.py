import os
from pathlib import Path
import pandas as pd
from namepipe import nt, NameTask

import kg_create_data
import kg_build_index

from kg_build_msa import (
    kirToMultiMsa,
    kirToSingleMsa,
    kirMerge2dl1s1,
)

from kg_utils import (
    runShell,
    runDocker,
    runAll,
    getSamples,
    threads,
    samtobam,
)

import kg_typping
import kg_typping_linnil1
from kg_eval import EvaluateKIR
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
def hisatMap(input_name, index):
    # 1 to 1
    output_name = input_name + "." + index.replace("/", "_")
    if Path(f"{output_name}.bam").exists():
        return output_name
    f1, f2 = input_name + ".read.1.fq", input_name + ".read.2.fq"
    runDocker("hisat", f"""\
              hisat2 --threads {threads} -x {index}.graph -1 {f1} -2 {f2} \
              --no-spliced-alignment --max-altstried 64 --haplotype \
              -S {output_name}.sam """)
    samtobam(output_name)
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
def hisatTyping(input_name, index):
    suffix = kg_typping.main(index, input_name + ".bam")
    return input_name + suffix


@nt
def hisatKIRTyping(input_name, index, exon=False, cohert=False):
    kir = kg_typping_linnil1.HisatKIR(index)
    if cohert:
        kir.cn_cohert = True

    if exon:
        kir.cn_type = "sam_exon_depth"
    else:
        kir.cn_type = "sam_depth"

    # kir.typing_by = "hisat"
    # kir.cn_dev = 0.04
    kir.cn_kde = False
    kir.cn_median = False
    # kir.typing_by = "hisat"
    # kir.typing_by = "likelihood_multi"
    kir.typing_by = "likelihood"
    suffix = kir.main(input_name)
    print(input_name, suffix)
    return input_name + suffix


@nt
def hisatKIRCohertCN(input_name, index):
    names = input_name.get_input_names()
    kir = kg_typping_linnil1.HisatKIR(index)
    # must be same as hisatKIRTyping
    # kir.cn_dev = 0.04
    kir.cn_kde = True
    kir.cn_median = True
    kir.cn_type = "sam_depth"     
    kir.typing_by = "likelihood"

    print(names)
    gene_cns = kir.calcAndSaveCohertCN(names, input_name.replace_wildcard("_mergeforCN"))
    print(gene_cns)

    return input_name.replace_wildcard("_mergeforCN")


@nt
def hisatKIRRessult(input_name, answer):
    ans = EvaluateKIR(f"{answer}/{answer}.summary.csv")

    called_alleles_dict = {}
    predict_list = []

    for name in input_name.get_input_names():
        predict = kg_typping_linnil1.readAlleleResult(name + ".tsv")[0]
        id = name.template_args[0]
        predict['id'] = id
        predict_list.append({
            'id': id,
            'name': predict['name'],
            'alleles': '_'.join(predict['alleles']),
        })
        called_alleles_dict[id] = predict['alleles']

    df = pd.DataFrame(predict_list)
    output_name = input_name.replace_wildcard("_merge")
    df.to_csv(f"{output_name}.tsv", index=False, sep="\t")
    print(df)
    ans.compareCohert(called_alleles_dict, skip_empty=True)
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


if __name__ == "__main__":
    data_folder = "data"
    Path(data_folder).mkdir(exist_ok=True)
    extract_exon = True

    answer_folder = "linnil1_syn_wide"
    answer_folder = "linnil1_syn_exon"
    answer_folder = "linnil1_syn_30x"
    answer_folder = "linnil1_syn_30x_seed87"
    Path(answer_folder).mkdir(exist_ok=True)

    index_folder = "index"
    os.makedirs(index_folder, exist_ok=True)

    samples = answer_folder >> createSamples

    if extract_exon:
        samples = samples >> extractExon.set_args(folder=answer_folder + "_exon")
        runShell(f"ln -fs ../{answer_folder}/{answer_folder}.summary.csv {answer_folder}_exon/{answer_folder}_exon.summary.csv")
        answer_folder = answer_folder + "_exon"

    samples = samples >> linkSamples.set_args(data_folder)
    if answer_folder == "linnil1_syn_wide":
        samples = samples >> link10Samples

    # msa_index = index_folder >> NameTask(func=kirToMultiMsa)  # "index/kir_2100_raw.mut01"
    msa_index = index_folder >> NameTask(func=kirMerge2dl1s1) # index = "index/kir_2100_2dl1s1.mut01"
    # msa_index = index_folder >> NameTask(func=kirToMultiMsa).set_args(split_2DL5=True) # index = "index/kir_2100_ab.mut01"
    # msa_index = index_folder >> NameTask(func=kirToSingleMsa) # index = "index/kir_2100_merge.mut01"

    index = msa_index >> NameTask(func=kg_build_index.main)
    print(index)
    # samples = "data2/linnil1_syn_30x_seed87.{}"
    mapping = samples >> hisatMap.set_args(index=str(index)) >> hisatTyping.set_args(index=str(index))
    typing = mapping >> hisatKIRTyping.set_args(index=str(index), exon=extract_exon)
    # typing_cohert = mapping >> hisatKIRCohertCN.set_args(index=str(index)).set_depended(0)
    # typing =        mapping >> hisatKIRTyping.set_args(index=str(index), exon=extract_exon, cohert=True)
    typing = typing >> hisatKIRRessult.set_args(answer=answer_folder).set_depended(0)
    # print(mapping)


    # bowtie mapping rate
    # bowtie2_index = index >> bowtie2BuildConsensus  # "index/kir_2100_?_cons"
    # bowtie2_index = index >> "index/kir_2100_raw.mut01" >> bowtie2BuildFull >> "index/kir_2100_raw_full" # index = "index/kir_2100_raw_full"
    # mapping = samples >> bowtie2.set_args(index=str(bowtie2_index))
    # mapping = samples >> bowtie2.set_args(index=str(bowtie2_index), use_arg="ping")

    # bwa
    # bwa_index = index >> bwaIndex  # "index/kir_2100_?_cons"
    # mapping = samples >> bwa.set_args(index=str(bwa_index))

    """
    index = "index/kir_2100_muscle"
    # kirToSingleMsa(method='muscle')
    # bamFilter("-f 0x2", ".nosingle")
    # suffix += ".nosec"
    # bamFilter("-f 0x2 -F 256", ".sec")
    """
    runShell("stty echo opost")
