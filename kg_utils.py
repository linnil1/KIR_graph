import os
import sys
import time
import uuid
import pickle
import traceback
import __main__
import subprocess
import importlib.util
from glob import glob
from typing import Callable, Any
from pathlib import Path

from Bio import SeqIO
from namepipe import BaseTaskExecutor, NamePath, StandaloneTaskExecutor
from graphkir.utils import (
    runShell,
    runDocker as runDockerGK,
)
from kg_eval import compareCohort, readPredictResult, readAnswerAllele


images = {
    'samtools': "quay.io/biocontainers/samtools:1.15.1--h1170115_0",
    'clustalo': "quay.io/biocontainers/clustalo:1.2.4--h1b792b2_4",
    'hisat':    "quay.io/biocontainers/hisat2:2.2.1--h87f3376_4",
    "muscle":   "quay.io/biocontainers/muscle:5.1--h9f5acd7_1",
    "bowtie":   "quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0",
    "picard":   "quay.io/biocontainers/picard:2.27.3--hdfd78af_0",
    "gatk3":    "docker.io/broadinstitute/gatk3:3.6-0",
    "gatk4":    "docker.io/broadinstitute/gatk:4.2.6.1",
    "bwa":      "quay.io/biocontainers/bwa:0.7.17--hed695b0_7",
    "deepvariant": "docker.io/google/deepvariant:1.4.0",
    "kpi":      "docker.io/droeatumn/kpi",
    "bcftools": "quay.io/biocontainers/bcftools:1.13--h3a49de5_0",
    "vg":       "quay.io/vgteam/vg:v1.40.0",
    "ping":     "localhost/linnil1/ping",
    "t1k":      "localhost/linnil1/t1k",
    "kpi":      "localhost/linnil1/kpi",
    "node":     "docker.io/library/node:18-alpine",
}


def runDocker(image: str, cmd: str, *arg, **kwargs):
    """ run docker container """
    image = images.get(image, image)
    return runDockerGK(image, cmd, *arg, **kwargs)


def buildDocker(image: str, dockerfile: str, folder: str = "."):
    return runShell(f"podman build {folder} -f {dockerfile} -t {images[image]}")


class SlurmTaskExecutor(StandaloneTaskExecutor):
    """ Run the task in SLURM """

    def __init__(
        self,
        template_file: str = "taiwania.template",
        threads_per_sample: int = 2,
        filename_template: str ="slurm_io/job_{time}_{input_name}",
    ):
        super().__init__(threads=1, filename_template=filename_template)
        self.template = open(template_file).read()
        self.threads_per_sample = threads_per_sample

    def submit_standalone_task(self, name: str) -> subprocess.CompletedProcess:
        """
        Save the command into shell script(xx.tmp.sh).

        And submit it to SLURM.

        The return process is the process to submit job (<1s),
        so, it will not block when task is not finished.
        """
        cmd = ["python << EOF",
               "from namepipe import StandaloneTaskExecutor",
               "from graphkir.utils import setThreads",
               f"setThreads({self.threads_per_sample})",
               "StandaloneTaskExecutor.run_standalone_task("
               f"{repr(str(__main__.__file__))}, {repr(str(name))})",
               f"EOF"]

        text = self.template.format(cmd="\n".join(cmd), name=name)
        with open(f"{name}.tmp.sh", "w") as f:
            f.write(text)
        return runShell(f"sbatch {name}.tmp.sh")

    def cleanup_task(self, name: str):
        """ remove script and log """
        super().cleanup_task(name)
        Path(f"{name}.tmp.sh").unlink()
        # Path(f"{name}.std.log").unlink()
        # Path(f"{name}.err.log").unlink()


def linkSamples(input_name, data_folder, new_name=None,
                fastq=True, fasta=False, sam=False):
    # 1 -> 1
    # only work for unix relative folder
    if new_name is None:
        name = Path(input_name.template).name
    else:
        name = new_name
    Path(data_folder).mkdir(exist_ok=True)
    output_name = str(Path(data_folder) / name)
    new_name = output_name.format(input_name.template_args[0])
    relative = "../" * (len(Path(new_name).parents) - 1)
    if fastq:
        if Path(new_name + ".read.1.fq").exists():
            return output_name
        runShell(f"ln -s {relative}{input_name}.read.1.fq {new_name}.read.1.fq")  # add .read is for previous naming
        runShell(f"ln -s {relative}{input_name}.read.2.fq {new_name}.read.2.fq")
    if fasta:
        if Path(new_name + ".fa").exists():
            return output_name
        runShell(f"ln -s {relative}{input_name}.fa {new_name}.fa")
    if sam:
        if Path(new_name + ".sam").exists():
            return output_name
        if Path(f"{input_name}.read..sam").exists():
            runShell(f"ln -s {relative}{input_name}.read..sam {new_name}.sam")
        else:
            runShell(f"ln -s {relative}{input_name}.sam {new_name}.sam")
    return output_name


def getAnswerFile(sample_name: str) -> str:
    """ The answer of cohort xxx.{}.oo is located at xxx_summary.oo.csv """
    # TODO: tempoary solution
    sample_name = NamePath(sample_name)
    name = sample_name.replace_wildcard("_summary")
    if Path(name + ".csv").exists():
        return name + ".csv"
    elif Path(name + ".tsv").exists():
        return name + ".tsv"
    else:
        raise ValueError(f"Not found answer file {name}")


def compareResult(input_name, sample_name, input_fasta_name=None, plot=False):
    print(input_name + ".tsv")
    answer_seq = {}
    # back = xx.00.30x.fq -> xx.00.fa
    for name in NamePath(back(NamePath(sample_name))).get_input_names():
        if not Path(name + ".fa").exists():
            continue
        # read answer allele
        # brute-force to rewrite "xxx-1" -> "xxx"
        allele_seq = SeqIO.to_dict(SeqIO.parse(name + ".fa", "fasta"))
        for i in list(allele_seq):
            if i.endswith("-1") or i.endswith("-2") or i.endswith("-3"):
                allele_seq[i[:-2]] = allele_seq[i]

        answer_seq[name.template_args[-1]] = allele_seq

    predit_seq = {}
    if input_fasta_name:
        for name in NamePath(input_fasta_name).get_input_names():
            if Path(name + ".fa").exists():
                print("HI", name + ".fa")
                predit_seq[name.template_args[-1]] = SeqIO.to_dict(SeqIO.parse(name + ".fa", "fasta"))

    compareCohort(
        readAnswerAllele(getAnswerFile(sample_name)),
        readPredictResult(input_name + ".tsv"),
        cohort_answer_seqs=answer_seq,
        cohort_predit_seqs=predit_seq,
        skip_empty=True,
        # base_compare=True,
        plot=plot,
    )
    return input_name


def addSuffix(input_name, suffix):
    return input_name + suffix


def back(input_name):
    return str(Path(input_name.template).with_suffix(""))
