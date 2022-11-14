"""
Utilities
* run shell
* run docker and image's version
* sam to bam
* numpy to json
"""
import json
import uuid
import subprocess
import dataclasses
import numpy as np


threads = 14
docker_path = "podman"
images = {
    'samtools': "quay.io/biocontainers/samtools:1.15.1--h1170115_0",
    'clustalo': "quay.io/biocontainers/clustalo:1.2.4--h1b792b2_4",
    'hisat':    "quay.io/biocontainers/hisat2:2.2.1--h87f3376_4",
    "muscle":   "quay.io/biocontainers/muscle:5.1--h9f5acd7_1",
}


def runDocker(image: str, cmd: str, capture_output=False, cwd=None, opts="") -> subprocess.CompletedProcess:
    """ run docker container """
    image = images.get(image, image)
    name = str(uuid.uuid4()).split("-", 1)[0]
    # docker_path = "docker"
    proc = runShell(f"{docker_path} run -it --rm -u root --name {name} "
                    f"-w /app -v $PWD:/app {opts} {image} {cmd}",
                    capture_output=capture_output, cwd=cwd)
    return proc


def runShell(cmd: str, capture_output=False, cwd=None) -> subprocess.CompletedProcess:
    """ wrap os.system """
    print(cmd)
    proc = subprocess.run(cmd, shell=True,
                          capture_output=capture_output,
                          cwd=cwd,
                          check=True,
                          universal_newlines=True)
    return proc


def samtobam(name: str, keep=False) -> None:
    """ samfile -> sorted bamfile and index (This is so useful) """
    runDocker("samtools", f"samtools sort -@{threads} {name}.sam -o {name}.bam")
    runDocker("samtools", f"samtools index            {name}.bam")
    if not keep:
        runShell(f"rm {name}.sam")


class NumpyEncoder(json.JSONEncoder):
    """ The encoder for saving numpy array to json """
    def default(self, obj):
        if dataclasses.is_dataclass(obj):
            return dataclasses.asdict(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def getGeneName(allele: str) -> str:
    """ KIR3DP1*BACKBONE -> KIR3DP1 """
    return allele.split("*")[0]


def limitAlleleField(allele: str, resolution: int = 7) -> str:
    """ KIR3DP1*0010101 with resolution 5 -> KIR3DP1*00101 """
    return getGeneName(allele) + "*" + getAlleleField(allele, resolution)


def getAlleleField(allele: str, resolution: int = 7) -> str:
    """ KIR3DP1*0010101 with resolution 5 -> 00101 """
    if "*" not in allele:
        return ""
    return allele.split("*")[1][:resolution]
