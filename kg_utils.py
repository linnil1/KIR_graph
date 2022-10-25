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

from namepipe import BaseTaskExecutor, NamePath
from graphkir.utils import (
    runShell,
    runDocker as runDockerGK,
)

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


class IsolateTaskExecutor(BaseTaskExecutor):
    """ Run the task in different process """

    def __init__(self, threads: int = 100):
        super().__init__()
        self.executor = ('python -c "from kg_utils import IsolateTaskExecutor; '
                         'IsolateTaskExecutor.runIsolateInCommandLine()" '
                         f'{__main__.__file__}')
        self.threads = threads

    @staticmethod
    def wildcardImportFromFile(path: str) -> None:
        """ same as from main.py import * """
        spec = importlib.util.spec_from_file_location("tmp", path)
        assert spec
        assert spec.loader
        module = importlib.util.module_from_spec(spec)
        assert module
        spec.loader.exec_module(module)
        for i in dir(module):
            if not i.startswith("__"):
                setattr(__main__, i, getattr(module, i))

    @staticmethod
    def runIsolateFunc(func: Callable, input_name: str) -> Any:
        """ run the function """
        try:
            return func(input_name)
        except BaseException as e:
            print(e)
            return BaseException(traceback.format_exc())

    @staticmethod
    def runIsolate(main_file_path: str, name: str) -> None:
        """ read input, run the function and save the output """
        IsolateTaskExecutor.wildcardImportFromFile(main_file_path)
        func, input_name = pickle.load(open(f"{name}.tmp.in", "rb"))
        out = IsolateTaskExecutor.runIsolateFunc(func, input_name)
        pickle.dump(out, open(f"{name}.tmp.out", "wb"), pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def runIsolateInCommandLine() -> None:
        """ same as runIsolateFunc but read sys.argv """
        if sys.argv[1] == "-c":
            main_file_path, name = sys.argv[2], sys.argv[3]
        else:
            main_file_path, name = sys.argv[1], sys.argv[2]
        IsolateTaskExecutor.runIsolate(main_file_path, name)

    def submitIsolateTask(self, name: str) -> subprocess.CompletedProcess:
        """ Define how to run the task in isolated python """
        cmd = f"{self.executor} {name}"
        return runShell(cmd)

    def submit(self, func: Callable, input_name: NamePath) -> str:
        """
        Setup the data needed for isolate task.
        * name
        * {name}.tmp.in
        * {name}.tmp.out
        """
        name = "job_" + str(input_name).replace("/", "_")
        self._logger.info(f"Write {name}.tmp.in")
        # create tmp input parameters
        pickle.dump(
            (func, input_name), open(f"{name}.tmp.in", "wb"), pickle.HIGHEST_PROTOCOL
        )
        # submit/start each task
        self.submitIsolateTask(name)
        return f"{name}.tmp.out"

    def run_tasks(self, names: list[NamePath], func: Callable) -> list[NamePath | str]:
        """
        Run the isolated tasks in threads pool to avoid blocking.

        If task preparion or task submittion fail, it'll immediatly fail.
        Otherwise, it will raise error after all processes are done.
        """
        # for name in names:
        #     output_tmp_files.append(self.submit(func, name))
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            exes = [executor.submit(self.submit, func, name) for name in names]
            self._logger.info(f"Submit {len(exes)} tasks")
            output_tmp_files = [exe.result() for exe in exes]
            self._logger.info(f"{len(exes)} tasks done")

        # use xx.tmp.out to check if the task is done
        output_name = []
        for output_tmp_file in output_tmp_files:
            while not os.path.exists(output_tmp_file):
                time.sleep(3)
            self._logger.info(f"Load {output_tmp_file}")
            output_name.append(pickle.load(open(output_tmp_file, "rb")))

        # print and raise error if return Error
        has_error = False
        for i in output_name:
            if isinstance(output_name[-1], BaseException):
                self._logger.error(f"{output_name[-1]}")
                has_error = True
        if has_error:
            raise ValueError("Some tasks have error")
        return output_name


class SlurmTaskExecutor(IsolateTaskExecutor):
    """ Run the task in SLURM """

    def __init__(self, template_file="taiwania.template"):
        super().__init__(threads=1)
        self.template = open(template_file).read()

    def submitIsolateTask(self, name: str) -> subprocess.CompletedProcess:
        """
        Save the command into shell script(xx.tmp.sh).

        And submit it to SLURM.

        The return process is the process to submit job (<1s),
        so, it will not block when task is not finished.
        """
        text = self.template.format(cmd=f"{self.executor} {name}", name=name)
        open(f"{name}.tmp.sh", "w").write(text)
        return runShell(f"sbatch {name}.tmp.sh")
