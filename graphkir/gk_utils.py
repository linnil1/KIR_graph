import os
from glob import glob
import uuid
import subprocess
from concurrent.futures import ProcessPoolExecutor


threads = 20
docker_path = "podman"
images = {
    'samtools': "quay.io/biocontainers/samtools:1.15.1--h1170115_0",
    'clustalo': "quay.io/biocontainers/clustalo:1.2.4--h1b792b2_4",
    'hisat':    "quay.io/biocontainers/hisat2:2.2.1--h87f3376_4",
    "muscle":   "quay.io/biocontainers/muscle:5.1--h9f5acd7_1",
}


def runDocker(image, cmd, capture_output=False):
    """ run docker container """
    if image in images:
        image = images[image]
    name = str(uuid.uuid4()).split("-")[0]
    # docker_path = "docker"
    proc = runShell(f"{docker_path} run -it --rm -u root --name {name} "
                    f"-w /app -v $PWD:/app {image} {cmd}",
                    capture_output=capture_output)
    return proc


def runShell(cmd, capture_output=False, cwd=None):
    """ wrap os.system """
    print(cmd)
    proc = subprocess.run(cmd, shell=True,
                          capture_output=capture_output,
                          cwd=cwd,
                          universal_newlines=True)
    proc.check_returncode()
    return proc


def samtobam(name, keep=False):
    """ This is so useful """
    runDocker("samtools", f"samtools sort -@4 {name}.sam -o {name}.bam")
    runDocker("samtools", f"samtools index    {name}.bam")
    if not keep:
        runShell(f"rm {name}.sam")
