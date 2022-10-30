from pathlib import Path
from functools import partial
import pandas as pd

from namepipe import compose, NameTask, ConcurrentTaskExecutor
from graphkir.utils import runShell, threads
from kg_utils import runDocker
from kg_eval import compareCohort, readPredictResult, readAnswerAllele
from kg_main import linkSamples, getAnswerFile


def downloadT1k(input_name):
    folder = "T1K"
    if Path(f"{folder}/genotyper").exists():
        return folder
    if not Path(folder).exists():
        runShell(f"git clone https://github.com/mourisl/T1K {folder}")
    runShell(f"docker build . -f t1k.dockerfile -t c4lab/t1k")
    runDocker("c4lab/t1k", f"make -j {threads}", cwd=folder)
    if not Path(f"{folder}/hlaidx").exists():
        runDocker("c4lab/t1k",
                  "perl t1k-build.pl -o hlaidx --download IPD-IMGT/HLA",
                  cwd=folder)
    if not Path(f"{folder}/kiridx").exists():
        runDocker("c4lab/t1k",
                  "perl t1k-build.pl -o kiridx --download IPD-KIR",
                  cwd=folder)
    return folder


def runT1k(input_name, index):
    output_name = input_name + ".t1k"
    if Path(output_name + "._genotype.tsv").exists():
        return output_name
    runDocker(
        "c4lab/t1k",
        f""" \
      {index}/run-t1k \
      -1 {input_name}.read.1.fq \
      -2 {input_name}.read.2.fq \
      --preset kir-wgs -f {index}/kiridx/kiridx_dna_seq.fa \
      -t {threads} \
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


def mergeT1kResult(input_name):
    output_name = input_name.replace_wildcard("_mergecall")
    called_alleles_dict = {}
    for name in input_name.list_names():
        alleles = readT1kAllele(name + "._genotype.tsv")
        id = name.template_args[0]
        called_alleles_dict[id] = alleles

    predict_list = []
    for id, alleles in called_alleles_dict.items():
        predict_list.append(
            {
                "id": id,
                "alleles": "_".join(alleles),
                "name": input_name.format(id),
            }
        )
    df = pd.DataFrame(predict_list)
    df.to_csv(f"{output_name}.tsv", index=False, sep="\t")
    print(df)
    return output_name


def compareT1kResult(input_name, sample_name):
    answer_file = getAnswerFile(sample_name)
    answer = readAnswerAllele(answer_file)
    predit = readPredictResult(input_name + ".tsv")
    print(input_name + ".tsv")
    compareCohort(answer, predit, skip_empty=True)
    return input_name


if __name__ == "__main__":
    NameTask.set_default_executor(ConcurrentTaskExecutor(threads=20))
    index_t1k = downloadT1k("")
    samples = "linnil1_syn/linnil1_syn_s44.{}.30x_s444"
    data_folder = "data5"

    compose([
            samples,
            partial(linkSamples, data_folder=data_folder),
            partial(runT1k, index=str(index_t1k)),
            NameTask(mergeT1kResult, depended_pos=[-1]),
            partial(compareT1kResult, sample_name=samples),
    ])
