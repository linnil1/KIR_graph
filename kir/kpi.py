# https://github.com/droeatumn/kpi
from typing import Any
from pathlib import Path

import pandas as pd

from .kir_pipe import KirPipe


class KPI(KirPipe):
    name = "kpi"

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)
        self.version = "v1.1.1"
        self.images = {
            # "kpi": f"localhost/c4lab/kpi:{self.version}",
            "kpi": "docker.io/droeatumn/kpi",
        }
        self.folder_name = "kpi"

    def download(self, folder_base: str = "") -> str:
        """Download KPI dockerfile and build it"""
        folder = self.folder_name + "_" + self.escapeName(self.version)
        if Path(folder).exists():
            return folder
        self.runShell(f"git clone https://github.com/droeatumn/kpi.git {folder}")
        # if not self.checkImage("kpi"):
        #     self.buildImage("kpi", f"{folder}/Dockerfile", folder=folder, args={"version": self.version})
        return folder

    def run(self, input_name: str, index: str) -> str:
        """Run KPI (It used next flow)"""
        mapping_file = self.replaceWildcard(input_name, "_kpidatalist")
        if Path(mapping_file + ".txt").exists():
            return input_name + ".kpi_prediction"

        with open(mapping_file + ".txt", "w") as f:
            for name in self.listFiles(input_name):
                f1, f2 = f"{name}.read.1.fq", f"{name}.read.2.fq"
                if not Path(f1).exists():
                    f1, f2 = f"{name}.read.1.fq.gz", f"{name}.read.2.fq.gz"
                basename = Path(name).name
                suffix = ".kpi_" + self.escapeName(index)
                print(basename + ".kpi", f1, sep="\t", file=f)
                print(basename + ".kpi", f2, sep="\t", file=f)

        folder = Path(input_name).parents[0]
        self.runDocker(
            "kpi", f"/opt/kpi/main.nf --map {mapping_file}.txt --output {folder}"
        )
        return input_name + ".kpi_prediction"

    def mergeResult(self, input_name: str, index: str) -> str:
        """Merge KPI haplotype result to allele and cn"""
        haps = pd.read_csv(f"{index}/input/haps.txt", sep="\t")
        output_name_cn = self.replaceWildcard(input_name, "_merge_cn")
        output_name = self.replaceWildcard(input_name, "_merge_guess_allele")

        # if Path(output_name + ".tsv").exists():
        #     return output_name

        cn = {}
        guess_allele = []
        for name in self.listFiles(input_name):
            df = pd.read_csv(f"{name}.txt", sep="\t")
            haplo = df["haplotypes"][0]
            print(name, haplo)
            haplo = haplo.split("|")[0]
            df = haps[haps["nomenclature"].isin(haplo.split("+"))]

            df = df.drop(
                columns=["haplotype", "nomenclature", "Jiang 2012 freq", "structure"]
            )
            df = df.set_axis(map(lambda i: f"KIR{i}", df.columns), axis=1)

            # cn
            id = self.getID(name)
            cn[id] = dict(df.sum(axis=0))

            # allele
            alleles = []
            for gene_name, gene_cn in cn[id].items():
                alleles.extend([gene_name] * int(gene_cn))
            guess_allele.append(
                {
                    "id": id,
                    "alleles": "_".join(alleles),
                    "name": name,
                }
            )
        assert cn

        # cn
        df_cn = pd.DataFrame(cn)
        df_cn = df_cn.reset_index().rename(columns={"index": "gene"})
        print(df_cn)
        df_cn.to_csv(output_name_cn + ".csv", index=False)

        # allele
        df_allele = pd.DataFrame(guess_allele)
        df_allele.to_csv(f"{output_name}.tsv", index=False, sep="\t")
        print(df_allele)
        return output_name

    def runAll(self, samples: str) -> str:
        """Run all the script(Don't use this when building pipeline"""
        index = self.download()
        samples = self.run(samples, index=index)
        samples = self.mergeResult(samples, index=index)
        return samples
