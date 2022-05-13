import os
from glob import glob
import uuid
import subprocess


threads = 30
docker_path = "podman"
images = {
    'samtools': "quay.io/biocontainers/samtools:1.15.1--h1170115_0",
    'clustalo': "quay.io/biocontainers/clustalo:1.2.4--h1b792b2_4",
    'hisat':    "quay.io/biocontainers/hisat2:2.2.1--h87f3376_4",
    "muscle":   "quay.io/biocontainers/muscle:5.1--h9f5acd7_1",
    "bowtie":   "quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0",
}


def runDocker(image, cmd, capture_output=False):
    """ run docker container """
    if image in images:
        image = images[image]
    name = str(uuid.uuid4()).split("-")[0]
    # docker_path = "docker"
    return runShell(f"{docker_path} run -it --rm -u root --name {name} "
                    f"-w /app -v $PWD:/app {image} {cmd}",
                    capture_output=capture_output)


def runShell(cmd, capture_output=False):
    """ wrap os.system """
    print(cmd)
    proc = subprocess.run(cmd, shell=True,
                          capture_output=capture_output,
                          universal_newlines=True)
    proc.check_returncode()
    return proc


def getSamples(name, endswith="", strict=True, return_name=False):
    """
    Args:
      name(str): `folder/xxx.a.b`
      endswith(str): `.c`
      strict(bool):  If strict, the format `folder/xxx.a.b.{xx}.c` where `xx` has no `.`
    """
    files = glob(name + ".*")
    # check name
    for i in files:
        assert i.startswith(name + ".")
    if endswith:
        files = [i for i in files if i.endswith(endswith)]
    if strict:
        files = [i for i in files if i[len(name) + 1: -len(endswith)].count(".") == 0]
    if not return_name:
        return sorted(files)

    names = [i[len(name) + 1: -len(endswith)] for i in files]
    return sorted(zip(files, names))