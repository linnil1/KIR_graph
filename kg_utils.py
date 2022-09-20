import os
from glob import glob
import uuid
import subprocess
from concurrent.futures import ProcessPoolExecutor

from graphkir.gk_utils import (
    runShell,
    samtobam,
    runDocker as runDockerGK,
    getGeneName,
    limitAlleleField,
    getAlleleField,
)

threads = 20
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
