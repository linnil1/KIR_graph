from typing import Any
from pathlib import Path

import pandas as pd

from .kir_pipe import KirPipe


class PING(KirPipe):
    name = "ping"

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)
        self.version = "20220527"
        self.images = {
            "ping": f"localhost/c4lab/ping:{self.version}",
        }
        self.folder_name = "PING"

    def download(self, folder_base: str = "") -> str:
        """Download index"""
        folder = self.folder_name + self.escapeName(self.version)
        if Path(folder).exists():
            return folder
        self.runShell(f"git clone https://github.com/wesleymarin/PING.git {folder}")
        self.runShell("git checkout 4cd8592", cwd=folder)
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
            f"Rscript kir/ping.run.R",
            opts=f""" \
                -v $PWD/{index}/Resources:/app/Resources:ro \
                -e RAW_FASTQ_DIR={folder_in} \
                -e FASTQ_PATTERN=fq \
                -e THREADS={self.getThreads()} \
                -e RESULTS_DIR={folder_out} \
            """,
        )
        return folder_out

    def migrateSample(self, input_name: str) -> str:
        """Copy file to new directory"""
        folder = self.replaceWildcard(input_name, "_pinglist")
        if Path(folder).exists():
            return folder
        self.runShell(f"mkdir -p {folder}")
        for name in self.listFiles(input_name):
            f1, f2 = f"{name}.read.1.fq", f"{name}.read.2.fq"
            if not Path(f1).exists():
                f1, f2 = f"{name}.read.1.fq.gz", f"{name}.read.2.fq.gz"
            self.runShell(f"ln -s ../../{f1} {folder}/{Path(f1).name}")
            self.runShell(f"ln -s ../../{f2} {folder}/{Path(f2).name}")
        return folder

    def mergeResult(self, input_name: str) -> str:
        output_name = input_name + ".merge"
        # if Path(output_name + ".tsv").exists():
        #     return output_name
        data = self.readAllele(f"{input_name}/finalAlleleCalls.csv")

        result = []
        for name, alleles in data.items():
            result.append(
                {
                    "id": name,
                    "alleles": "_".join(alleles),
                    "name": name,
                }
            )
        assert result

        df = pd.DataFrame(result)
        df.to_csv(f"{output_name}.tsv", index=False, sep="\t")
        print(df)
        return output_name

    def readAllele(self, csv_file: str) -> dict[str, list[str]]:
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
        data = data.rename(columns={"Unnamed: 0": "name"})

        called_alleles = {}
        for sample in data.to_dict("records"):
            id = sample["name"]
            alleles = []
            for gene, alleles_str in sample.items():
                if gene == "name":
                    continue
                alleles.extend(alleles_str.split(" ")[0].split("+"))
            alleles = [i for i in alleles if "null" not in i]
            alleles = [i for i in alleles if "failed" not in i]
            # alleles = [i for i in alleles if "unresolved" not in i]
            called_alleles[id] = alleles
            # print(id, alleles)
        # print(called_alleles)
        return called_alleles

    def runAll(self, samples: str) -> str:
        """Run all the script(Don't use this when building pipeline)"""
        index = self.download()
        samples = "example_data/test.{}"
        samples = self.migrateSample(samples)
        # If you break at this line. Fine
        # Try again after editing manualCopyNumberFrame.csv
        samples = self.main(samples, index=index)
        samples = self.mergeResult(samples)
        return samples

    @classmethod
    def readGeneDepthRatio(cls, locus_csv: str) -> pd.DataFrame:
        """Read PING's gene depth format (locusRatioFrame)"""
        ping_data = pd.read_csv(locus_csv)
        ping_data = ping_data.rename(columns={"Unnamed: 0": "sample"})
        ping_data["method"] = "PING"
        ping_data["id"] = list(ping_data["sample"])
        ping_data = ping_data.drop(columns=["sample"])
        return ping_data
