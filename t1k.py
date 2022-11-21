from pathlib import Path
from functools import partial
from itertools import chain
import pandas as pd

from namepipe import compose, NameTask, ConcurrentTaskExecutor
from graphkir.utils import runShell, getThreads
from kg_utils import runDocker, linkSamples, getAnswerFile, compareResult
from kg_eval import saveCohortAllele


def downloadT1k(input_name):
    folder = "T1K"
    if Path(f"{folder}/genotyper").exists():
        return folder
    if not Path(folder).exists():
        runShell(f"git clone https://github.com/mourisl/T1K {folder}")
    runShell(f"docker build . -f t1k.dockerfile -t c4lab/t1k")
    runDocker("c4lab/t1k", f"make -j {getThreads()}", cwd=folder)
    if not Path(f"{folder}/hlaidx").exists():
        runDocker("c4lab/t1k",
                  "perl t1k-build.pl -o hlaidx --download IPD-IMGT/HLA",
                  cwd=folder)
    if not Path(f"{folder}/kiridx").exists():
        runDocker("c4lab/t1k",
                  "perl t1k-build.pl -o kiridx --download IPD-KIR",
                  cwd=folder)
    return folder


def runT1k(input_name, index, digits="7"):
    output_name = input_name + f".t1k_{digits}"
    if Path(output_name + "._genotype.tsv").exists():
        return output_name
    # relax = --relaxIntronAlign
    runDocker(
        "c4lab/t1k",
        f""" \
      {index}/run-t1k \
      -1 {input_name}.read.1.fq \
      -2 {input_name}.read.2.fq \
      --preset kir-wgs -f {index}/kiridx/kiridx_dna_seq.fa \
      --alleleDigitUnits {digits} \
      -t {getThreads()} \
      -o {output_name}. \
    """,
    )
    return output_name


def readT1kAllele(t1k_tsv: str) -> list[str]:
    column = ["gene_name", "num_alleles",
              "allele_1", "abundance_1", "quality_1",
              "allele_2", "abundance_2", "quality_2"]
    df = pd.read_csv(t1k_tsv, sep="\t", names=column)
    # reorganize
    df1 = df[["allele_1", "abundance_1", "quality_1"]]
    df1.columns = ["allele", "abundance", "quality"]  # type: ignore
    df2 = df[["allele_2", "abundance_2", "quality_2"]]
    df2.columns = ["allele", "abundance", "quality"]  # type: ignore
    df = pd.concat([df1, df2])

    # remove low quality
    df = df[df["quality"] > 5]
    print(df)
    """
    allele  abundance  quality
    0    KIR2DL1*052  66.463100       23
    1    KIR2DL2*003  74.674877       27
    3    KIR2DL4*054  77.410274       58
    4   KIR2DL5A*001  47.363073       18
    """
    return list(df["allele"])


def mergeT1kResult(input_name, select="first"):
    output_name = input_name.replace_wildcard("_mergecall")
    if select == "all":
        output_name += "_all"
    elif select == "first":
        pass
    else:
        raise ValueError

    called_alleles_dict = {}
    for name in input_name.list_names():
        alleles = readT1kAllele(name + "._genotype.tsv")
        if select == "first":
            alleles = [i.split(",")[0] for i in alleles]
        elif select == "all":
            alleles = list(chain.from_iterable(i.split(",") for i in alleles))
        id = name.template_args[0]
        called_alleles_dict[id] = alleles

    saveCohortAllele(called_alleles_dict, f"{output_name}.tsv")
    return output_name


if __name__ == "__main__":
    NameTask.set_default_executor(ConcurrentTaskExecutor(threads=20))
    index_t1k = downloadT1k("")
    samples = "linnil1_syn/linnil1_syn_s44.{}.30x_s444"
    data_folder = "data6"
    # samples = "linnil1_syn/linnil1_syn_s2022.{}.30x_s1031"
    # data_folder = "data6"
    compose([
        samples,
        partial(linkSamples, data_folder=data_folder),
        partial(runT1k, index=str(index_t1k), digits=7),
        NameTask(mergeT1kResult, func_kwargs=dict(select="all"), depended_pos=[-1]),
        partial(compareResult, sample_name=samples),
    ])
    runShell("stty echo opost")
