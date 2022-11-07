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
from concurrent.futures import ThreadPoolExecutor

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
}


def runDocker(image: str, cmd: str, *arg, **kwargs):
    """ run docker container """
    image = images.get(image, image)
    return runDockerGK(image, cmd, *arg, **kwargs)


def runSingularity(
    image: str, cmd: str, capture_output=False, cwd=None
) -> subprocess.CompletedProcess:
    """ run siguarity which image is downloaded as file """
    image = images.get(image, image)
    return runShell(
        "singularity run "
        "--bind /staging/biology "
        "--bind /staging/biodata/lions/twbk/ "
        "--bind /work/linnil1tw/ "
        f"{image.split('/')[-1]}.sif "
        f"{cmd}",
        capture_output=capture_output,
    )


class SlurmTaskExecutor(StandaloneTaskExecutor):
    """ Run the task in SLURM """

    def __init__(self, template_file="taiwania.template"):
        super().__init__(threads=1)
        self.template = open(template_file).read()

    def submit_standalone_task(self, name: str) -> subprocess.CompletedProcess:
        """
        Save the command into shell script(xx.tmp.sh).

        And submit it to SLURM.

        The return process is the process to submit job (<1s),
        so, it will not block when task is not finished.
        """
        cmd = ["python", "<<", "EOF\n",
               "from namepipe import StandaloneTaskExecutor\n",
               "StandaloneTaskExecutor.run_standalone_task("
               f"{repr(str(__main__.__file__))}, {repr(str(name))})\n",
               f"EOF\n"]

        text = self.template.format(cmd=" ".join(cmd), name=name)
        with open(f"{name}.tmp.sh", "w") as f:
            f.write(text)
        return runShell(f"sbatch {name}.tmp.sh")


def linkSamples(input_name, data_folder, new_name=None, fastq=True, fasta=False, sam=False):
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
    return sample_name.replace_wildcard("_summary") + ".csv"


def compareResult(input_name, sample_name):
    print(input_name + ".tsv")
    print(readPredictResult(input_name + ".tsv"))
    compareCohort(
        readAnswerAllele(getAnswerFile(sample_name)),
        readPredictResult(input_name + ".tsv"),
        skip_empty=True,
        # plot=True,
    )
    return input_name


def addSuffix(input_name, suffix):
    return input_name + suffix


def back(input_name):
    return str(Path(input_name.template).with_suffix(""))
