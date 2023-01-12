"""T1K pipeline"""
from typing import Any
from pathlib import Path
from itertools import chain

import pandas as pd

from .kir_pipe import KirPipe


class T1k(KirPipe):
    """T1K pipeline https://github.com/mourisl/T1K"""

    name = "t1k"

    def __init__(self, version: str = "v1.0.1", **kwargs: Any):
        super().__init__(**kwargs)
        self.version = version
        self.images = {
            "t1k": f"localhost/c4lab/t1k:{version}",
        }
        self.setIPDVersion("2100")

    def download(self, folder_base: str = "") -> str:
        """Download t1k and compile"""
        if folder_base and not folder_base.endswith("/"):
            folder_base += "/"
        folder = folder_base + "t1k_" + self.escapeName(self.version)
        if Path(folder).exists():
            return folder
        if not self.checkImage("t1k"):
            self.buildImage("t1k", "kir/t1k.dockerfile", args={"version": self.version})
        self.runShell(f"mkdir -p {folder}")
        return folder

    def build(self, folder: str, ipd_version: str = "2100") -> str:
        """Build KIR database"""
        output_name = f"{folder}/ipd_{ipd_version}"
        if Path(output_name).exists():
            return output_name
        self.runShell(f"mkdir -p {output_name}")
        if ipd_version <= "2100":
            self.runShell(
                f"wget https://github.com/ANHIG/IPDKIR/raw/{ipd_version}/KIR.dat"
                f"  -O {output_name}/kir.dat"
            )
        else:
            self.runShell(
                f"wget https://github.com/ANHIG/IPDKIR/raw/{ipd_version}/kir.dat "
                f"  -O {output_name}/kir.dat"
            )
        self.runDocker(
            "t1k",
            f"t1k-build.pl -o {output_name} -d {output_name}/kir.dat --prefix t1k",
        )
        return output_name

    def run(self, input_name: str, index: str, digits: int = 7) -> str:
        """Main function"""
        output_name = input_name + f".t1k_{self.escapeName(index)}.dig{digits}"
        if Path(output_name + "._genotype.tsv").exists():
            return output_name
        # relax = --relaxIntronAlign
        f1, f2 = f"{input_name}.read.1.fq", f"{input_name}.read.2.fq"
        if not Path(f1).exists():
            f1, f2 = f"{input_name}.read.1.fq.gz", f"{input_name}.read.2.fq.gz"
        self.runDocker(
            "t1k",
            f""" \
            run-t1k \
              -1 {f1} -2 {f2} \
              --preset kir-wgs -f {index}/t1k_dna_seq.fa \
              --alleleDigitUnits {digits} \
              -t {self.getThreads()} \
              -o {output_name}. \
            """,
        )
        return output_name

    def readAlleles(self, t1k_tsv: str) -> list[str]:
        """Read T1k csv"""
        column = [
            "gene_name", "num_alleles",
            "allele_1", "abundance_1", "quality_1",
            "allele_2", "abundance_2", "quality_2",
        ]
        df = pd.read_csv(t1k_tsv, sep="\t", names=column)
        # reorganize
        df1 = df[          ["allele_1", "abundance_1", "quality_1"]]
        df1 = df1.set_axis(["allele",   "abundance",   "quality"], axis=1)
        df2 = df[          ["allele_2", "abundance_2", "quality_2"]]
        df2 = df2.set_axis(["allele",   "abundance",   "quality"], axis=1)
        df = pd.concat([df1, df2])

        # remove low quality
        df = df[df["quality"] > 5]
        # print(df)
        # allele  abundance  quality
        # 0    KIR2DL1*052  66.463100       23
        # 1    KIR2DL2*003  74.674877       27
        # 3    KIR2DL4*054  77.410274       58
        # 4   KIR2DL5A*001  47.363073       18
        return list(df["allele"])

    def mergeResult(self, input_name: str, select: str = "first") -> str:
        """Merge all samples result"""
        output_name = self.replaceWildcard(input_name, "_mergecall")
        if select == "all":
            output_name += "_all"
        elif select == "first":
            pass
        else:
            raise ValueError
        # if Path(output_name + ".tsv").exists():
        #     return output_name

        predict_list = []
        for name in self.listFiles(input_name):
            alleles = self.readAlleles(name + "._genotype.tsv")
            if select == "first":
                alleles = [i.split(",")[0] for i in alleles]
            elif select == "all":
                alleles = list(chain.from_iterable(i.split(",") for i in alleles))
            predict_list.append(
                {
                    "id": self.getID(name),
                    "alleles": alleles,
                    "name": name,
                }
            )
        self.savePredictedAllele(predict_list, output_name)
        return output_name

    def runAll(self, input_name: str) -> str:
        """Run all the script(Don't use this when building pipeline"""
        index = self.download()
        index = self.build(index, self.ipd_version)
        samples = input_name
        for sample in self.listFiles(samples):
            sample = self.run(sample, index=index)
        sample = self.mergeResult(samples + f".t1k_{self.escapeName(index)}.dig7")
        return sample
