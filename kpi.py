# https://github.com/droeatumn/kpi
import os
from pathlib import Path
import pandas as pd
from namepipe import nt, NameTask
from kg_utils import runDocker, runShell


@nt
def setup(input_name, kpi_folder):
    if Path(kpi_folder).exists():
        return kpi_folder
    runShell(f"git clone https://github.com/droeatumn/kpi.git {kpi_folder}")
    # runShell("docker build . -t kirkpi ", cwd=kpi_folder)
    return kpi_folder


@nt
def linkSamples(input_name, data_folder):
    # 1 -> 1
    name = input_name.split('/')[0]
    output_name = os.path.join(data_folder, name + ".{}")
    new_name = output_name.format(input_name.template_args[0])
    if Path(new_name + ".read.1.fq").exists():
        return output_name
    # using hard link because we change data folder
    runShell(f"ln {input_name}.1.fq {new_name}.read.1.fq")
    runShell(f"ln {input_name}.2.fq {new_name}.read.2.fq")
    return output_name


@nt
def runKPI(input_name, kpi_folder):
    mapping_file = input_name.replace_wildcard("_merge_for_mapping")
    if Path(mapping_file + ".txt").exists():
        return input_name + ".kpi_prediction"

    with open(mapping_file + ".txt", "w") as f:
        for name in input_name.get_input_names():
            basename = Path(name).name
            print(basename + ".kpi", name + ".read.1.fq",sep="\t", file=f)
            print(basename + ".kpi", name + ".read.2.fq",sep="\t", file=f)
    folder = Path(input_name).parents[0]

    runDocker("kpi", f"./main.nf --map {mapping_file}.txt --output {folder}",
              opts=f"-v $PWD/{folder}:/opt/kpi/{folder}", workdir="/opt/kpi")
    return input_name + ".kpi_prediction"


@nt
def collectResult(input_name, kpi_folder):
    haps = pd.read_csv(f"{kpi_folder}/input/haps.txt", sep="\t")

    output_name = input_name.replace_wildcard("_merge_cn")

    cn = []
    for name in input_name.get_input_names():
        df = pd.read_csv(f"{name}.txt", sep="\t")
        haplo = df['haplotypes'][0]
        print(name, haplo)
        haplo = haplo.split("|")[0]
        a = haps[haps["nomenclature"].isin(haplo.split("+"))]
        a = a.drop(columns=["haplotype", "nomenclature", "Jiang 2012 freq", "structure"])
        a.columns = map(lambda i: "KIR" + i, a.columns)
        gene_cn = dict(a.sum(axis=0))
        gene_cn['name'] = name
        cn.append(gene_cn)

    cn = pd.DataFrame(cn)
    cn = cn.set_index("name").astype(int)
    print(cn)
    cn.to_csv(output_name + '.csv')
    return output_name


if __name__ == "__main__":
    kpi_folder = "kirkpi"
    kpi = None >> setup.set_args(kpi_folder)
    answer = "linnil1_syn_30x"
    answer = "linnil1_syn_30x_seed87"
    data_folder = "data3"
    Path(data_folder).mkdir(exist_ok=True)
    samples = answer + "/" + answer + ".{}.read" >> linkSamples.set_args(data_folder)
    samples >> runKPI.set_depended(-1).set_args(kpi_folder) >> collectResult.set_depended(-1).set_args(kpi_folder)
    # ignore pseduo gene in kg_eval.py
