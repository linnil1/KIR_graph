import os
import json
from pathlib import Path
import pandas as pd
from pyhlamsa import KIRmsa, Genemsa

from graphkir.utils import runShell, samtobam, getThreads
from kg_utils import runDocker


def filterMSAname(input_name):
    return "KIR" in input_name.template_args[-1]


def runVG(cmd):
    return runDocker("vg", f"bash -c '{cmd}'")


def msa2vcf(input_name):
    if not filterMSAname(input_name):
        return input_name
    if Path(f"{input_name}.vcf.gz").exists():
        return input_name
    msa = Genemsa.load_msa(f"{input_name}.fa", f"{input_name}.json")
    msa.to_fasta(f"{input_name}.backbone.fa", gap=False, ref_only=True)
    msa.to_vcf(f"{input_name}.vcf", plain_text=True)
    runDocker("bcftools", f"bcftools sort {input_name}.vcf -o {input_name}.vcf.gz")
    runDocker("bcftools", f"bcftools index -t {input_name}.vcf.gz")
    return input_name


def buildVG(input_name):
    output_name = input_name.replace_wildcard("_vgindex")
    if Path(f"{output_name}.giraffe.gbz").exists():
        return output_name
    msa_names = [name for name in input_name.get_input_names() if filterMSAname(name)]
    runVG(f"""\
          vg autoindex -w giraffe \
          {"".join(f"-r {name}.backbone.fa " for name in msa_names)} \
          {"".join(f"--vcf {name}.vcf.gz " for name in msa_names)} \
          -p {output_name} """)
    return output_name


def mapVG(input_name, index, output_type="bam"):
    assert output_type in ["bam", "gam"]
    output_name = input_name + "." + index.replace("/", "_")
    if output_type == "bam" and Path(f"{output_name}.bam").exists():
        return output_name
    if output_type == "gam" and Path(f"{output_name}.gam").exists():
        return output_name
    runVG(f""" \
          vg giraffe -Z {index}.giraffe.gbz \
                     -m {index}.min -d {index}.dist \
                     -t {getThreads()} --fragment-mean 400 --fragment-stdev 2 \
                     -f {input_name}.read.1.fq -f {input_name}.read.2.fq \
          """ + (
                   f" -o SAM > {output_name}.sam        " if output_type == "bam" else
                   f" -o gam > {output_name}.unsort.gam "
          )
    )
    if output_type == "bam":
        samtobam(output_name)
    else:
        runVG(f"vg gamsort -t {getThreads()} {output_name}.unsort.gam "
              f"-i {output_name}.gam.gai > {output_name}.gam")
    return output_name


def vgindex2VGTube(input_name):
    output_name = input_name + ".vis"
    if Path(f"{output_name}.xg").exists():
        return output_name
    runVG(f"vg convert {input_name}.giraffe.gbz -x > {output_name}.xg")
    return output_name


def setupVGTube(input_name):
    if os.path.exists("sequenceTubeMap"):
        return "sequenceTubeMap"
    runShell("git clone https://github.com/vgteam/sequenceTubeMap.git")
    runShell("wget https://github.com/vgteam/vg/releases/download/v1.42.0/vg")
    runShell("chmod +x vg")
    runShell("mv vg sequenceTubeMap")
    runDocker("node", "yarn",         cwd="sequenceTubeMap")
    runDocker("node", "yarn install", cwd="sequenceTubeMap")
    runDocker("node", "yarn build",   cwd="sequenceTubeMap")
    return "sequenceTubeMap"


def writeVGTubeSettings(input_name, index, tube_path):
    index_folder = os.path.split(index)[0]
    data_folder = os.path.split(input_name)[0]
    runShell(f"ln -fs ../{index_folder} {tube_path}/{index_folder}")
    runShell(f"ln -fs ../{data_folder} {tube_path}/{data_folder}")
    data = {
        "BACKEND_URL": False,
        "DATA_SOURCES": [
            {
                "name": name,
                "xgFile": f"{index}.xg",
                "gamFile": f"{name}.gam",
                "dataPath": "mounted",
                "region": "KIR3DL3*BACKBONE:1000-2000",
                "dataType": "built-in",
            } for name in input_name.get_input_names()],
        "vgPath": "",
        "dataPath": ".",
        "internalDataPath": "",
    }
    # from pprint import pprint
    # pprint(data)
    json.dump(data, open(f"{tube_path}/src/config.json", "w"))
    runDocker("node", "yarn build", cwd=tube_path)
    return input_name


def startVGTube(input_name, tube_path, port=8000):
    runShell(f"docker run -it --rm -p {port}:3000 "
             f"-v $PWD:/app -w /app/{tube_path} -v $PWD/{tube_path}/vg:/bin/vg  "
             "node:18-alpine yarn start")
    return tube_path


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
