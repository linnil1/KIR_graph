# https://github.com/droeatumn/kpi
from pathlib import Path
from functools import partial
import pandas as pd

from namepipe import nt, NamePath
from graphkir.utils import runShell
from kg_eval import saveCohortAllele
from kg_utils import runDocker, buildDocker, linkSamples, compareResult


def setup(kpi_folder):
    if Path(kpi_folder).exists():
        return kpi_folder
    runShell(f"git clone https://github.com/droeatumn/kpi.git {kpi_folder}")
    buildDocker("kpi", f"{kpi_folder}/Dockerfile", folder=kpi_folder)
    return kpi_folder


@nt(depended_pos=[-1])
def runKPI(input_name, kpi_folder):
    mapping_file = input_name.replace_wildcard("_merge_for_mapping")
    if Path(mapping_file + ".txt").exists():
        return input_name + ".kpi_prediction"

    with open(mapping_file + ".txt", "w") as f:
        for name in input_name.get_input_names():
            basename = Path(name).name
            print(basename + ".kpi", name + ".read.1.fq", sep="\t", file=f)
            print(basename + ".kpi", name + ".read.2.fq", sep="\t", file=f)

    folder = Path(input_name).parents[0]
    runDocker("kpi", f"/opt/kpi/main.nf --map {mapping_file}.txt --output {folder}")
    return input_name + ".kpi_prediction"


@nt(depended_pos=[-1])
def collectResult(input_name, kpi_folder):
    haps = pd.read_csv(f"{kpi_folder}/input/haps.txt", sep="\t")
    output_name_cn = input_name.replace_wildcard("_merge_cn")
    output_name = input_name.replace_wildcard("_merge_guess_allele")

    cn = []
    guess_allele = {}
    for name in input_name.get_input_names():
        df = pd.read_csv(f"{name}.txt", sep="\t")
        haplo = df['haplotypes'][0]
        print(name, haplo)
        haplo = haplo.split("|")[0]
        a = haps[haps["nomenclature"].isin(haplo.split("+"))]
        a = a.drop(columns=["haplotype", "nomenclature", "Jiang 2012 freq", "structure"])
        a.columns = map(lambda i: "KIR" + i, a.columns)

        # cn
        gene_cn = dict(a.sum(axis=0))
        gene_cn['name'] = name
        cn.append(gene_cn)
        print(gene_cn)

        # allele
        guess_allele[name.template_args[0]] = []
        for gene_name, gene_cn in gene_cn.items():
            if gene_name == "name":
                continue
            for _ in range(int(gene_cn)):
                guess_allele[name.template_args[0]].append(gene_name)

    # cn
    cn = pd.DataFrame(cn)
    cn = cn.set_index("name").astype(int)
    print(cn)
    cn.to_csv(output_name_cn + '.csv')

    # allele
    saveCohortAllele(guess_allele, output_name + ".tsv")
    return output_name


if __name__ == "__main__":
    samples = "linnil1_syn/linnil1_syn_s44.{}.30x_s444"
    data_folder = "data6"
    kpi_folder = "kirkpi"
    setup(kpi_folder=kpi_folder)
    NamePath(samples) \
            >> partial(linkSamples, data_folder=data_folder) \
            >> runKPI(func_kwargs=dict(kpi_folder=kpi_folder)) \
            >> collectResult(func_kwargs=dict(kpi_folder=kpi_folder)) \
            >> partial(compareResult, sample_name=samples)
    # ignore pseduo gene in kg_eval.py
