"""
Utilities
* run shell
* run docker and image's version
* sam to bam
* numpy to json
"""
from glob import glob
from typing import TypedDict, Any
import re
import json
import uuid
import logging
import subprocess
import dataclasses

import numpy as np
import pandas as pd
from pyhlamsa import Genemsa

import requests
from tqdm import tqdm


class _ResType(TypedDict):
    threads: int
    memory: int
    engine: str


# logging
# What the fuck python logging is still hard to used
# logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("graphkir")
logger.propagate = False
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(logging.Formatter("%(asctime)s [%(name)s] [%(levelname)8s] %(message)s"))
logger.addHandler(ch)


def downloadFile(url: str, output_path: str) -> None:
    """Download file using requests with progress bar"""
    logger.info(f"Downloading {url} to {output_path}")
    response = requests.get(url, stream=True)
    response.raise_for_status()

    total_size = int(response.headers.get("content-length", 0))
    block_size = 8192

    with open(output_path, "wb") as f:
        with tqdm(
            total=total_size, unit="B", unit_scale=True, desc="Downloading"
        ) as pbar:
            for chunk in response.iter_content(chunk_size=block_size):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))


resources: _ResType = {  # per sample
    "threads": 2,
    "memory": 7,  # unit: G
    "engine": "podman",
}
engine_config = {
    "podman": {
        "path": "podman",
        "run": "run -t --rm -u root -w /app -v $PWD:/app",
        "name": "--name",  # short name has collision when submit the jobs concurrently
        "image_func": lambda i: i,
    },
    "docker": {
        "path": "/usr/bin/docker",
        "run": "run -t --rm -u root -w /app -v $PWD:/app",
        "name": "--name",
        "image_func": lambda i: i,
    },
    "singularity": {
        "path": "singularity",
        "run": "run -B $PWD ",
        "image_func": lambda i: "docker://" + i,
    },
    "local": {
        "path": "",
        "run": "",
        "image_func": lambda i: "",
    },
    # linnil1's configuration in Taiwania HPC
    "singularity_linnil1": {
        "path": "singularity",
        "run": (
            "run "
            "--bind /home/linnil1tw "
            "--bind /work/linnil1tw "
            "--bind /staging/biology/linnil1tw "
            "--bind /staging/biology/zxc898977 "
            "--bind /staging/biodata/lions/twbk "
        ),
        # lions: twbb data
        # zxc898977: hprc data
        "image_func": lambda i: "singur_image/" + i,
    },
}

# download every image into image file
# singularity build quay.io/biocontainers/bwa:0.7.17--hed695b0_7 \
#      singur_image/quay.io/biocontainers/bwa:0.7.17--hed695b0_7
images = {
    "samtools": "quay.io/biocontainers/samtools:1.15.1--h1170115_0",
    "clustalo": "quay.io/biocontainers/clustalo:1.2.4--h1b792b2_4",
    "hisat":    "quay.io/biocontainers/hisat2:2.2.1--h87f3376_4",
    "muscle":   "quay.io/biocontainers/muscle:5.1--h9f5acd7_1",
    "bwa":      "quay.io/biocontainers/bwa:0.7.17--hed695b0_7",
}


def getThreads() -> int:
    """Get number of thread per job"""
    return resources["threads"]


def setThreads(threads: int) -> None:
    """Set number of thread per job"""
    global resources
    resources["threads"] = threads


def getEngine() -> str:
    """Get engine(podman/docker/local/singularity)"""
    return resources["engine"]


def setEngine(engine: str) -> None:
    """Set engine(podman/docker/local/singularity)"""
    global resources
    resources["engine"] = engine
    assert engine in engine_config


def runDocker(
    image: str,
    cmd: str,
    capture_output: bool = False,
    cwd: str | None = None,
    opts: str = "",
) -> subprocess.CompletedProcess[str]:
    """run docker container"""
    image = images.get(image, image)
    random_name = str(uuid.uuid4()).split("-", 1)[0]

    # bad but works
    if "singularity" in getEngine():
        opts = opts.replace(" -e ", " --env ")
        opts = opts.replace(" -v ", " --bind ")
        opts = opts.replace(":/app/", ":$PWD/")
    conf = engine_config[getEngine()]
    cmd_all = (
        f"{conf['path']} {conf['run']} "
        + (f"{conf['name']} {random_name} " if conf.get("name") else "")
        + f"{opts} {conf['image_func'](image)} {cmd}"  # type: ignore
    )
    proc = runShell(cmd_all, capture_output=capture_output, cwd=cwd)
    return proc


def runShell(
    cmd: str, capture_output: bool = False, cwd: str | None = None
) -> subprocess.CompletedProcess[str]:
    """wrap os.system"""
    logger.debug(f"[Run] {cmd}")
    proc = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        cwd=cwd,
        check=True,
        universal_newlines=True,
    )
    if not capture_output:
        logger.debug(proc.stdout)
    return proc


def samtobam(name: str, keep: bool = False) -> None:
    """samfile -> sorted bamfile and index (This is so useful)"""
    runDocker("samtools", f"samtools sort -@{getThreads()}  {name}.sam -o {name}.bam")
    runDocker("samtools", f"samtools index -@{getThreads()} {name}.bam")
    if not keep:
        runShell(f"rm {name}.sam")


class NumpyEncoder(json.JSONEncoder):
    """The encoder for saving numpy array to json"""

    def default(self, obj: Any) -> Any:  # I don't know the exact format
        if dataclasses.is_dataclass(obj):
            return dataclasses.asdict(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def getGeneName(allele: str) -> str:
    """KIR3DP1*BACKBONE -> KIR3DP1"""
    return allele.split("*")[0]


def limitAlleleField(allele: str, resolution: int = 7) -> str:
    """KIR3DP1*0010101 with resolution 5 -> KIR3DP1*00101"""
    return getGeneName(allele) + "*" + getAlleleField(allele, resolution)


def getAlleleField(allele: str, resolution: int = 7) -> str:
    """
    KIR3DP1*0010101 with resolution 5 -> 00101
    KIR2DL1*0320102N with resolution 5 -> 003201
    KIR2DL1*0320102N with resolution 7 -> 00320102N
    """
    if "*" not in allele:
        return ""
    patterns = re.findall(r"^\w+\*(\d+\w*)", allele)
    if patterns:
        num = str(patterns[0])
    else:
        num = "new"
    if resolution == 7:
        return num
    else:
        return num[:resolution]
    # mostly equvient to this
    # return allele.split("*")[1][0][:resolution]


def mergeAllele(allele_result_files: list[str], final_result_file: str) -> pd.DataFrame:
    """Merge allele calling result"""
    df = pd.concat(pd.read_csv(f, sep="\t") for f in allele_result_files)
    df.to_csv(final_result_file, index=False, sep="\t")
    return df


def mergeCN(cn_result_files: list[str], final_result_file: str) -> pd.DataFrame:
    """Merge copy number result"""
    dfs = []
    for f in cn_result_files:
        df = pd.read_csv(f, sep="\t")
        df["name"] = f
        dfs.append(df)
    df = pd.pivot_table(pd.concat(dfs), values="cn", index="gene", columns=["name"])
    df = df.fillna(0)
    df = df.astype(int)
    df.to_csv(final_result_file, sep="\t")
    return df


def readFromMSAs(prefix: str) -> dict[str, Genemsa]:
    """
    Read MSAs via `{prefix}.*.json`

    Returns:
        genes: A dictionary of gene's name and its MSA
    """
    genes = {}
    for filename in glob(prefix + ".*.json"):
        split_name = filename[len(prefix) + 1 :].split(".")
        # prefix.anotherprefix.*.json will not included
        if len(split_name) != 2:
            continue
        logger.debug(f"[MSA] Read {filename}")
        gene = split_name[0]
        genes[gene] = Genemsa.load_msa(filename[:-5] + ".fa", filename)
    return genes
