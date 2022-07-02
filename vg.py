import os
from pathlib import Path
import pandas as pd
from namepipe import nt, NameTask
from kg_utils import runDocker, runShell, threads, samtobam
from pyhlamsa import KIRmsa, msaio


@nt
def linkSamples(input_name, data_folder):
    # 1 -> 1
    name = input_name.split('/')[0]
    output_name = os.path.join(data_folder, name + ".{}")
    new_name = output_name.format(input_name.template_args[0])
    if Path(new_name + ".read.1.fq").exists():
        return output_name
    # using hard link because we change data folder
    runShell(f"ln -s ../{input_name}.1.fq {new_name}.read.1.fq")
    runShell(f"ln -s ../{input_name}.2.fq {new_name}.read.2.fq")
    return output_name


@nt
def linkMSA(input_name, folder_name):
    # run after any  (in main)
    # NameTask(func=kirToMultiMsa)
    # NameTask(func=kirMerge2dl1s1)
    # NameTask(func=kirToMultiMsa).set_args(split_2DL5=True)
    # NameTask(func=kirToSingleMsa)

    name = Path(input_name).name
    output_name = f"{folder_name}/{name}"
    output_template = f"{folder_name}/{Path(input_name.template).name}"
    if input_name.template_args[0] in ["400bp_match", "variation"] :
        return output_template
    if Path(f"{output_name}.json").exists():
        return output_template

    runShell(f"cp {input_name}.json {output_name}.json")
    runShell(f"cp {input_name}.fa   {output_name}.fa")
    return output_template


@nt
def msa2vcf(input_name):
    if Path(f"{input_name}.vcf.gz").exists():
        return input_name
    msa = msaio.load_msa(f"{input_name}.fa", f"{input_name}.json")
    msaio.to_vcf(msa, f"{input_name}.vcf", plain_text=True)
    msaio.to_fasta(msa, f"{input_name}.ref.fa", gap=False, ref_only=True)
    runDocker("bcftools", f"bcftools sort {input_name}.vcf -o {input_name}.vcf.gz")
    runDocker("bcftools", f"bcftools index -t {input_name}.vcf.gz")
    return input_name 


@nt
def buildVG(input_name):
    output_name = input_name.replace_wildcard("_vgindex").replace("save_vgindex", "vgindex")
    if Path(f"{output_name}.giraffe.gbz").exists():
        return output_name
    
    runShell(f"cat "
             + " ".join(f"{name}.ref.fa" for name in input_name.get_input_names())
             + f" > {output_name}.ref.fa")
    runDocker("vg", f"""\
             vg autoindex -w giraffe \
             -r {output_name}.ref.fa \
             {"".join(f"--vcf {name}.vcf.gz " for name in input_name.get_input_names())} \
             -p {output_name} """)
    return output_name


@nt
def mapVG(input_name, index):
    output_name = input_name + "." + index.replace("/", "_")
    runDocker("vg", f""" \
        vg giraffe -Z {index}.giraffe.gbz \
                   -m {index}.min -d {index}.dist \
                   -t {threads} --fragment-mean 400 --fragment-stdev 2 \
                   -f {input_name}.read.1.fq -f {input_name}.read.2.fq \
                   -o SAM > {output_name}.sam
    """)
    samtobam(output_name)
    return output_name


data_folder = "data3"
Path(data_folder).mkdir(exist_ok=True)
index = "index/kir_2100_raw.save.{}" >> linkMSA.set_args(data_folder)
index = "index/kir_2100_merge.save.{}" >> linkMSA.set_args(data_folder)
index = index >> msa2vcf >> buildVG.set_depended(-1)
print(index)

answer = "linnil1_syn_30x"
answer = "linnil1_syn_30x_seed87"
samples = answer + "/" + answer + ".{}.read" >> linkSamples.set_args(data_folder)
samples >> mapVG.set_args(index=str(index))
# ignore pseduo gene in kg_eval.py
