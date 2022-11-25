"""
Utilities
* run shell
* run docker and image's version
* sam to bam
* numpy to json
"""
import re
import json
import uuid
import subprocess
import dataclasses
import numpy as np

resources: dict[str, int] = {  # per sample
    'threads': 2,
    'memory': 7,  # unit: G
}
docker_config = {
    "podman": {
        "path": "podman",
        "run": "run -it --rm -u root -w /app -v $PWD:/app",
        "name": "--name",  # short name has collision when submit the jobs concurrently
        "image_prefix": "",
    },
    "docker": {
        "path": "/usr/bin/docker",
        "run": "run -it --rm -u root -w /app -v $PWD:/app",
        "name": "--name",
        "image_prefix": "",
    },
    "singularity": {
        "path": "singularity",
        "run": "run -B $PWD ",
        "image_prefix": "docker://",
    },
}
docker_type = "podman"
# linnil1's configuration in Taiwania HPC
# download every image into image file
# singularity build quay.io/biocontainers/bwa:0.7.17--hed695b0_7 singur_image/quay.io/biocontainers/bwa:0.7.17--hed695b0_7
# docker_type = "singularity"
# docker_config[docker_type]["image_prefix"] = "singur_image/"
# docker_config[docker_type]["run"] = (
#     "run "
#     "--bind /home/linnil1tw "
#     "--bind /work/linnil1tw "
#     "--bind /staging/biology/linnil1tw "
#     "--bind /staging/biology/zxc898977 "
#     "--bind /staging/biodata/lions/twbk "
# )
# lions: twbb data
# zxc898977: hprc data 
images = {
    "samtools": "quay.io/biocontainers/samtools:1.15.1--h1170115_0",
    "clustalo": "quay.io/biocontainers/clustalo:1.2.4--h1b792b2_4",
    "hisat":    "quay.io/biocontainers/hisat2:2.2.1--h87f3376_4",
    "muscle":   "quay.io/biocontainers/muscle:5.1--h9f5acd7_1",
}


def getThreads() -> int:
    return resources['threads']


def setThreads(threads: int) -> None:
    global resources
    resources['threads'] = threads


def runDocker(image: str, cmd: str, capture_output=False, cwd=None, opts="") -> subprocess.CompletedProcess:
    """ run docker container """
    image = images.get(image, image)
    random_name = str(uuid.uuid4()).split("-", 1)[0]

    conf = docker_config[docker_type]
    cmd_all = f"{conf['path']} {conf['run']} " + \
              (f"{conf['name']} {random_name} " if conf.get('name') else "") + \
              f"{opts} {conf['image_prefix']}{image} {cmd}"
    proc = runShell(cmd_all, capture_output=capture_output, cwd=cwd)
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
    runDocker("samtools", f"samtools sort -@{getThreads()}  {name}.sam -o {name}.bam")
    runDocker("samtools", f"samtools index -@{getThreads()} {name}.bam")
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
    """
    KIR3DP1*0010101 with resolution 5 -> 00101
    KIR2DL1*0320102N with resolution 5 -> 003201
    KIR2DL1*0320102N with resolution 7 -> 00320102N
    """
    if "*" not in allele:
        return ""
    patterns = re.findall(r"^\w+\*(\d+\w*)", allele)
    if patterns:
        num = patterns[0]
    else:
        num = "new"
    if resolution == 7:
        return num
    else:
        return num[:resolution]
    # return allele.split("*")[1][0][:resolution]
