import os
import json
from pathlib import Path
import pandas as pd
from namepipe import nt, NameTask
from pyhlamsa import KIRmsa, msaio

from kg_main import linkSamples
from kg_utils import runDocker, runShell, threads, samtobam


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
    if input_name.template_args[0] in ["400bp_match", "gene_variation", "variation"] :
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
    if Path(f"{output_name}.bam").exists():
        return output_name
    runDocker("vg", f""" \
        vg giraffe -Z {index}.giraffe.gbz \
                   -m {index}.min -d {index}.dist \
                   -t {threads} --fragment-mean 400 --fragment-stdev 2 \
                   -f {input_name}.read.1.fq -f {input_name}.read.2.fq \
                   -o SAM > {output_name}.sam
    """)
    samtobam(output_name)
    return output_name


@nt
def visualizeMapping(input_name, index):
    output_name = input_name + "." + index.replace("/", "_")
    if Path(f"{output_name}.gam").exists():
        return output_name
    runDocker("vg", f""" bash -c '\
        vg giraffe -Z {index}.giraffe.gbz \
                   -m {index}.min -d {index}.dist \
                   -t {threads} --fragment-mean 400 --fragment-stdev 2 \
                   -f {input_name}.read.1.fq -f {input_name}.read.2.fq \
                   -o gam > {output_name}.unsort.gam '""")
    runDocker("vg", f"bash -c 'vg gamsort -t {threads} {output_name}.unsort.gam "
                    f"-i {output_name}.gam.gai > {output_name}.gam'")
    return output_name


@nt
def visualizeIndex(input_name):
    output_name = input_name + ".vis"
    if Path(f"{output_name}.xg").exists():
        return output_name
    runDocker("vg", f"bash -c 'vg convert {input_name}.giraffe.gbz -x > {output_name}.xg'")
    return output_name


@nt
def setupTube(input_name):
    if os.path.exists("sequenceTubeMap"):
        return "sequenceTubeMap"
    runShell("git clone https://github.com/vgteam/sequenceTubeMap.git")
    runShell("wget https://github.com/vgteam/vg/releases/download/v1.42.0/vg")
    runShell("chmod +x vg")
    runShell("mv vg sequenceTubeMap")
    runDocker("node:18-alpine", "yarn",         cwd="sequenceTubeMap")
    runDocker("node:18-alpine", "yarn install", cwd="sequenceTubeMap")
    runDocker("node:18-alpine", "yarn build",   cwd="sequenceTubeMap")
    return "sequenceTubeMap"


@nt
def writeVisualizeSettings(input_name, index, tube_path):
    folder = os.path.split(input_name)[0]
    data = {
        "BACKEND_URL": False,
        "DATA_SOURCES": [
            {
                "name": name,
                "xgFile": f"{index}.xg",
                "gamFile": f"{name}.gam",
                "dataPath": "mounted",
                "region": "KIR3DL3*BACKBONE:1-1000",
                "dataType": "built-in",
            } for name in input_name.get_input_names()],
        "vgPath": "",
        "dataPath": ".",
        "internalDataPath": "",
    }
    # from pprint import pprint
    # pprint(data)
    json.dump(data, open(f"{tube_path}/src/config.json", "w"))
    runDocker("node:18-alpine", "yarn build", cwd=tube_path)
    runShell("docker run -it --rm -w /app -p 8000:3000 "
             f"-v $PWD/{tube_path}:/app -v $PWD/{tube_path}/vg:/bin/vg -v $PWD/{folder}:/app/{folder} "
             "node:18-alpine yarn start")

@nt
def index2Svg(input_name):
    # runShell("docker build . -f graphviz.dockerfile -t graphviz")
    # output_name = f"{input_name}.kir3dl3"
    # runDocker("vg", f"bash -c \"vg find -x {input_name}.xg -p 'KIR3DL3*BACKBONE' | vg view -dp - > {output_name}.dot\"")
    # output_name = f"{input_name}.all"
    # runDocker("vg", f"bash -c \"vg view -x {input_name}.xg -dp - > {output_name}.dot\"")
    # output_name = f"{input_name}.kir3dl3-1-10000"
    # runDocker("vg", f"bash -c \"vg find -x {input_name}.xg -p 'KIR3DL3*BACKBONE:1-10000 -c 100' | vg view -dp - > {output_name}.dot\"")
    # output_name = f"{input_name}.kir3dl1-kir3dl2"
    # runDocker("vg", f"bash -c \"vg find -x {input_name}.xg -p 'KIR3DL3*BACKBONE' -p 'KIR3DL2*BACKBONE' | vg view -dp - > {output_name}.dot\"")
    output_name = f"{input_name}.kir2dl"
    runDocker("vg", f"bash -c \"vg find -x {input_name}.xg -p 'KIR2DL1*BACKBONE' -p 'KIR2DL2*BACKBONE' -p 'KIR2DL3*BACKBONE' -p 'KIR2DL4*BACKBONE'  -p 'KIR2DL5*BACKBONE'| vg view -dp - > {output_name}.dot\"")
    runDocker("graphviz", f"dot -Tsvg {output_name}.dot -o {output_name}.svg")
    print(f"{output_name}.svg")


if __name__ == "__main__":
    data_folder = "data3"
    Path(data_folder).mkdir(exist_ok=True)
    # index = "index/kir_2100_merge.save.{}" >> linkMSA.set_args(data_folder)
    index = "index/kir_2100_raw.save.{}" >> linkMSA.set_args(data_folder)
    index = index >> msa2vcf >> buildVG.set_depended(-1)
    print(index)

    answer = "linnil1_syn_30x"
    answer = "linnil1_syn_30x_seed87"
    samples = answer + "/" + answer + ".{}.read" >> nt(linkSamples).set_args(data_folder)
    samples >> mapVG.set_args(index=str(index))

    # visualize
    index_vis = index >> visualizeIndex
    samples_vis = samples >> visualizeMapping.set_args(index=str(index))
    tube = "" >> setupTube
    # samples_vis >> writeVisualizeSettings.set_depended(-1).set_args(index=str(index_vis), tube_path=str(tube))
    index_vis >> index2Svg
