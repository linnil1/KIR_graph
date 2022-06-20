import os
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
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
from kg_ping_threshold_plot import cutThresholdByAns


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
    assert "kir_2100_raw" in name
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
def hisatTyping(input_name, index):
    suffix = kg_typping.main(index, input_name + ".bam")
    return input_name + suffix


@nt
def hisatKIRTyping(input_name, index, exon=False):
    kir = kg_typping_linnil1.HisatKIR(index)
    if exon:
        kir.cn_type = "sam_exon_depth"
        kir.cn_median = True
    else:
        kir.cn_type = "sam_depth"
    kir.typing_by = "likelihood"
    suffix = kir.main(input_name)
    print(input_name, suffix)
    return input_name + suffix


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
def buildPing(input_name):
    folder = "PING"
    if os.path.exists(folder):
        return folder
    runShell(f"git clone https://github.com/wesleymarin/PING.git {folder}")
    runShell(f"docker build {folder} -t ping")
    runShell(f"ln -s ../{folder} {folder}")
    return folder


@nt
def pingCopyFile(input_name, index="PING"):
    # generate isolated fastq folder
    name = Path(input_name).name
    folder_in = os.path.join(index, input_name.replace("/", "_").replace(f".{input_name.template_args[0]}", ""))
    Path(folder_in).mkdir(exist_ok=True)
    if Path(f"{folder_in}/{name}.read.1.fq").exists():
        return folder_in

    # ln -fs ../data2/linnil1_syn_30x_seed87.00.read.1.fq PING/data2_linnil1_syn_30x_seed87/linnil1_syn_30x_seed87.00.read.1.fq
    relative = "../" * len(Path(folder_in).parents)
    runShell(f"ln -s {relative}{input_name}.read.1.fq {folder_in}/{name}.read.1.fq")
    runShell(f"ln -s {relative}{input_name}.read.2.fq {folder_in}/{name}.read.2.fq")
    return folder_in


@nt
def ping(input_name, index):
    assert index == "PING"
    folder_in = input_name
    folder_out = input_name + ".result"

    # first time will fail
    # but we'll get locusRatioFrame.csv
    if not os.path.exists(folder_out + "/locusRatioFrame.csv"):
        try:
            # -e SHORTNAME_DELIM={threads}
            runDocker("ping", f"Rscript PING_run.R", opts=f""" \
                -e RAW_FASTQ_DIR={folder_in} \
                -e FASTQ_PATTERN=fq \
                -e THREADS={threads} \
                -e RESULTS_DIR={folder_out} \
            """)
        except subprocess.CalledProcessError:
            pass

        # Use answer and ratio to cut thresholds
        assert os.path.exists(folder_out + "/locusRatioFrame.csv")
        # data2_linnil1_syn_30x_seed87/
        # -> linnil1_syn_30x_seed87
        name = input_name.split('/')[-1].split('.')[0].split("_", maxsplit=1)[1]
        print(f"Set CN threshold fro PING", name)
        cutThresholdByAns(
            f"{name}/{name}.summary.csv",
            folder_out
        )

    # This will success
    if not os.path.exists(folder_out + "/finalAlleleCalls.csv"):
        runDocker("ping", f"Rscript PING_run.R", opts=f""" \
            -e RAW_FASTQ_DIR={folder_in} \
            -e FASTQ_PATTERN=fq \
            -e THREADS={threads} \
            -e RESULTS_DIR={folder_out} \
        """)
    return folder_out + "/finalAlleleCalls"


@nt
def pingResult(input_name, answer_name):
    kir = EvaluateKIR(f"{answer_name}/{answer_name}.summary.csv")
    ping_called_alleles = readPingResult(f"{input_name}.csv")
    kir.compareCohert(ping_called_alleles)
    return input_name


@nt
def extractExon(input_name, folder):
    Path(folder).mkdir(exist_ok=True)
    copy_name = f"{folder}/{Path(input_name).name.replace('.read', '')}"
    output_name  = copy_name + "_exon"

    if Path(f"{output_name}.read.1.fq").exists():
        return output_name.replace(input_name.template_args[0], "{}") + ".read"

    runShell(f"ln -s ../{input_name.replace('.read', '')}.fa {copy_name}.fa")
    runShell(f"ln -s ../{input_name}..sam                    {copy_name}.sam")
    extractExonPairReads(f"{copy_name}.fa", f"{copy_name}.sam", f"{output_name}")
    return output_name.replace(input_name.template_args[0], "{}") + ".read"


if __name__ == "__main__":
    data_folder = "data"
    Path(data_folder).mkdir(exist_ok=True)
    extract_exon = True

    answer_folder = "linnil1_syn_wide"
    answer_folder = "linnil1_syn_exon"
    answer_folder = "linnil1_syn_30x_seed87"
    answer_folder = "linnil1_syn_30x"
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

    msa_index = index_folder >> NameTask(func=kirToMultiMsa)  # "index/kir_2100_raw.mut01"
    # msa_index = index_folder >> NameTask(func=kirMerge2dl1s1) # index = "index/kir_2100_2dl1s1.mut01"
    # msa_index = index_folder >> NameTask(func=kirToMultiMsa).set_args(split_2DL5=True) # index = "index/kir_2100_ab.mut01"
    # msa_index = index_folder >> NameTask(func=kirToSingleMsa) # index = "index/kir_2100_merge.mut01"

    index = msa_index >> NameTask(func=kg_build_index.main)
    print(index)
    # samples = "data2/linnil1_syn_30x_seed87.{}"
    mapping = samples >> hisatMap.set_args(index=str(index)) >> hisatTyping.set_args(index=str(index)) >> hisatKIRTyping.set_args(index=str(index), exon=extractExon)
    mapping >> hisatKIRRessult.set_args(answer=answer_folder).set_depended(0)
    # print(mapping)

    # ping_index = None >> buildPing
    # ping_predict = samples >> pingCopyFile.set_args(index=str(ping_index)) >> ping.set_args(index=str(ping_index)) 
    # ping_predict >> pingResult.set_args(answer_name=answer_folder)
    # print(ping_predict)

    # bowtie mapping rate
    # bowtie2_index = index >> bowtie2BuildConsensus >> "index2/kir_2100_raw_cons"
    # bowtie2_index = index >> bowtie2BuildFull >> "index2/kir_2100_raw_full" # index = "index/kir_2100_raw_full"
    # bowtie2_index = index >> "index2/kir_2100_ab.mut01" >> bowtie2BuildConsensus >> "index2/kir_2100_ab_cons"
    # samples >> bowtie2.set_args(index=str(bowtie2_index))
    # samples >> bowtie2.set_args(index=str(bowtie2_index), use_arg="ping")
    """
    index = "index/kir_2100_muscle"
    # kirToSingleMsa(method='muscle')
    # bamFilter("-f 0x2", ".nosingle")
    # suffix += ".nosec"
    # bamFilter("-f 0x2 -F 256", ".sec")
    """
    runShell("stty echo opost")
