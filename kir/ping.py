"""PING pipeline (targeted and wgs version included)"""
from typing import Any
from pathlib import Path

import pandas as pd

from .kir_pipe import KirPipe, logger


class PING(KirPipe):
    """
    PING pipeline(targeted and wgs version included)
    https://github.com/wesleymarin/PING.git
    """

    name = "ping"

    def __init__(self, version: str = "20220527", **kwargs: Any):
        super().__init__(**kwargs)
        self.version = version
        self.images = {
            "ping": f"localhost/c4lab/ping:{self.version}",
        }

    def download(self, folder_base: str = "") -> str:
        """Download index"""
        if folder_base and not folder_base.endswith("/"):
            folder_base += "/"
        folder = folder_base + "ping_" + self.escapeName(self.version)
        if Path(folder).exists():
            return folder
        self.runShell(f"git clone https://github.com/wesleymarin/PING.git {folder}")
        if self.version == "20220527":
            self.runShell("git checkout 4cd8592", cwd=folder)
        elif self.version == "wgs":
            self.runShell("git checkout wgs_snakemake", cwd=folder)
        else:
            raise ValueError("Ping version not found")
        if not self.checkImage("ping"):
            self.buildImage("ping", "kir/ping.dockerfile")
        return folder

    def getOutputFolder(self, folder_in: str, index: str) -> str:
        """A handly function"""
        return folder_in + ".result_" + self.escapeName(index)

    def main(self, folder_in: str, index: str) -> str:
        """Run PING_run"""
        folder_out = self.getOutputFolder(folder_in, index)
        if Path(folder_out + "/finalAlleleCalls.csv").exists():
            return folder_out

        self.runDocker(
            "ping",
            f"Rscript kir/ping.run_{self.version}.R",
            opts=f""" \
                -v $PWD/{index}/Resources:/app/Resources:ro \
                -e RAW_FASTQ_DIR={folder_in} \
                -e FASTQ_PATTERN=fq \
                -e THREADS={self.getThreads()} \
                -e RESULTS_DIR={folder_out} \
                -e SHORTNAME_DELIM=.read \
            """,
        )
        return folder_out

    def migrateSample(self, input_name: str) -> str:
        """Copy file to new directory"""
        folder = self.replaceWildcard(input_name, "_pingsample")
        if Path(folder).exists():
            return folder
        self.runShell(f"mkdir -p {folder}")
        for name in self.listFiles(input_name):
            f1, f2 = f"{name}.read.1.fq", f"{name}.read.2.fq"
            suffix = "fq"
            if not Path(f1).exists():
                f1, f2 = f"{name}.read.1.fq.gz", f"{name}.read.2.fq.gz"
                suffix = "fq.gz"
            self.runShell(f"ln -s ../../{f1} {folder}/id.{self.getID(name)}.read.1.{suffix}")
            self.runShell(f"ln -s ../../{f2} {folder}/id.{self.getID(name)}.read.2.{suffix}")
        return folder

    def mergeResult(self, input_name: str, use_novel: bool = False) -> str:
        """Read PING result finalAlleleCalls"""
        output_name = input_name + ".merge"
        # if Path(output_name + ".tsv").exists():
        #     return output_name
        if use_novel:
            output_name += "_iter"
            data = self.readAllele(f"{input_name}/iterAlleleCalls.csv")
        else:
            output_name += "_final"
            data = self.readAllele(f"{input_name}/finalAlleleCalls.csv")
        logger.debug(f"[PING] Raw result {data}")

        predict_list = []
        for name, alleles in data.items():
            predict_list.append(
                {
                    "id": name,
                    "alleles": alleles,
                    "name": input_name + "." + name,
                }
            )
        self.savePredictedAllele(predict_list, output_name)
        return output_name

    @classmethod
    def readAllele(cls, csv_file: str) -> dict[str, list[str]]:
        """
        Read ping result finalAlleleCalls.csv

        Format:
        ```
        name                      KIR3DP1                                        KIR2DS35
        linnil1_syn_wide.00.read. KIR3DP1*026+KIR3DP1*null                       KIR2DS3*009+KIR2DS3*024+KIR2DS5*02701
        linnil1_syn_wide.01.read. KIR3DP1*00302+KIR3DP1*03201 KIR3DP1*00304+KIR3
        """
        # read
        data = pd.read_csv(csv_file, sep=",")
        if not isinstance(data.index, pd.RangeIndex):
            data = data.reset_index()
        data = data.rename(columns={"Unnamed: 0": "name", "index": "name"})
        data = data.fillna("")

        called_alleles = {}
        for sample in data.to_dict("records"):
            name_id = sample["name"][3:]  # remove id.
            alleles = []
            for gene, alleles_str in sample.items():
                if gene == "name":
                    continue
                alleles.extend(alleles_str.split(" ")[0].split("+"))
            alleles = [i for i in alleles if i]
            alleles = [i for i in alleles if "null" not in i]
            alleles = [i for i in alleles if "failed" not in i]
            alleles = [i.replace("_", ".") for i in alleles]
            # alleles = [i for i in alleles if "unresolved" not in i]
            called_alleles[name_id] = alleles
            # print(name_id, alleles)
        # print(called_alleles)
        return called_alleles

    def runAll(self, input_name: str) -> str:
        """Run all the script(Don't use this when building pipeline)"""
        logger.info("[PING] Download references")
        index = self.download()
        samples = input_name
        samples = self.migrateSample(samples)
        # If you break at this line. Fine
        # Try again after editing manualCopyNumberFrame.csv
        logger.info(f"[PING] Run {samples}")
        samples = self.main(samples, index=index)
        samples = self.mergeResult(samples)
        return samples

    @classmethod
    def readGeneDepthRatio(cls, locus_csv: str) -> pd.DataFrame:
        """Read PING's gene depth format (locusRatioFrame)"""
        ping_data = pd.read_csv(locus_csv)
        ping_data = ping_data.rename(columns={"Unnamed: 0": "sample"})
        ping_data["method"] = "PING"
        ping_data["id"] = list(
            map(lambda i: str(i)[3:], ping_data["sample"])
        )  # remove id.
        ping_data = ping_data.drop(columns=["sample"])
        return ping_data
