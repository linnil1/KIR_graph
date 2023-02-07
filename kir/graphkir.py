"""Graph-KIR (my command line)"""
from typing import Any
from pathlib import Path

import pandas as pd

from .kir_pipe import KirPipe, logger


class GraphKir(KirPipe):
    """Run command graphkir"""

    name = "graphkir"

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)
        self.version = "alpha"
        self.run_engine = self.executor.engine
        # self.run_engine = "singularity_linnil1"

    def download(self, folder_base: str = "") -> str:
        """Graphkir Build stage (Build from database)"""
        if folder_base and not folder_base.endswith("/"):
            folder_base += "/"
        folder = folder_base + "graphkir_" + self.escapeName(self.version)
        if Path(folder).exists():
            return folder
        self.runShell(
            f""" \
            graphkir \
                --thread {self.getThreads()} \
                --r1 data/linnil1_syn_s44.00.30x_s444.read.1.fq \
                --r2 data/linnil1_syn_s44.00.30x_s444.read.2.fq \
                --index-folder {folder} \
                --output-folder {folder}/tmp \
                --output-cohort-name {folder}/tmp/example_data.cohort \
                --step-skip-typing \
                --log-level DEBUG \
                --engine {self.run_engine}
            """
        )
        return folder

    def run(self, input_name: str, index: str) -> str:
        """Run graphkir command"""
        # samples to csv
        mapping_file = self.replaceWildcard(input_name, "_graphkir_list")
        out_suffix = ".graphkir_" + self.escapeName(index)
        output_name  = self.replaceWildcard(input_name, out_suffix.replace(".", "_") + "_predict_cohort")
        if Path(output_name + ".allele1.tsv").exists():
            return output_name
        prefix_to_id = {}
        with open(mapping_file + ".csv", "w") as f:
            print("name,r1,r2,id", file=f)
            for name in self.listFiles(input_name):
                f1, f2 = f"{name}.read.1.fq", f"{name}.read.2.fq"
                if not Path(f1).exists():
                    f1, f2 = f"{name}.read.1.fq.gz", f"{name}.read.2.fq.gz"
                if not Path(f1).exists():
                    print(f"[name].read.1.fq or [name].read.1.fq.gz NOT FOUND. Skip {name=}")
                    continue
                sample_id = self.getID(name)
                print(name + out_suffix, f1, f2, sample_id, sep=",", file=f)
                prefix_to_id[name + out_suffix] = sample_id

        # Input csv and output csv
        self.runShell(
            f""" \
            graphkir \
                --thread {self.getThreads()} \
                --input-csv {mapping_file}.csv \
                --index-folder {index} \
                --cn-individually \
                --allele-method exonfirst \
                --output-cohort-name {output_name} \
                --log-level DEBUG \
                --engine {self.run_engine}
            """
        )

        # Brute force find ID
        df = pd.read_csv(f"{output_name}.allele.tsv", sep="\t")
        names = df["name"]
        for name in names:
            for prefix, sample_id in prefix_to_id.items():
                if name.startswith(prefix):
                    df.loc[df["name"] == name, "id"] = sample_id
        df.to_csv(f"{output_name}.allele1.tsv", index=False, sep="\t")
        return output_name

    def runAll(self, input_name: str) -> str:
        """Run all the script"""
        logger.info("[GraphKir] Setup Index")
        index = self.download()
        samples = input_name
        logger.info(f"[GraphKir] Run {samples}")
        samples = self.run(samples, index=index)
        return samples + ".allele1"
