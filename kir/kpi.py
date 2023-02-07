"""KIR-KPI pipeline"""
from typing import Any
from pathlib import Path

import pandas as pd

from .kir_pipe import KirPipe, logger


class KPI(KirPipe):
    """KIR-KPI pipeline (https://github.com/droeatumn/kpi)"""

    name = "kpi"

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)
        self.version = "v1.1.1"
        self.images = {
            # "kpi": f"localhost/c4lab/kpi:{self.version}",
            "kpi": "docker.io/droeatumn/kpi",
        }

    def download(self, folder_base: str = "") -> str:
        """Download KPI dockerfile and build it"""
        if folder_base and not folder_base.endswith("/"):
            folder_base += "/"
        folder = folder_base + "kpi_" + self.escapeName(self.version)
        if Path(folder).exists():
            return folder
        self.runShell(f"git clone https://github.com/droeatumn/kpi.git {folder}")
        # if not self.checkImage("kpi"):
        #     self.buildImage("kpi", f"{folder}/Dockerfile", folder=folder, args={"version": self.version})
        return folder

    def run(self, input_name: str, index: str) -> str:
        """Run KPI (It used next flow)"""
        mapping_file = self.replaceWildcard(input_name, "_kpidatalist")
        out_suffix = ".kpi_" + self.escapeName(index)
        output_name = input_name + out_suffix + "_prediction"
        if Path(mapping_file + ".txt").exists():
            return output_name
        with open(mapping_file + ".txt", "w") as f:
            for name in self.listFiles(input_name):
                f1, f2 = f"{name}.read.1.fq", f"{name}.read.2.fq"
                if not Path(f1).exists():
                    f1, f2 = f"{name}.read.1.fq.gz", f"{name}.read.2.fq.gz"
                basename = Path(name).name
                print(basename + out_suffix, f1, sep="\t", file=f)
                print(basename + out_suffix, f2, sep="\t", file=f)

        folder = Path(input_name).parents[0]
        self.runDocker(
            "kpi", f"/opt/kpi/main.nf --map {mapping_file}.txt --output {folder}"
        )
        return output_name

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
            logger.debug(f"[KPI] {name} {haplo}")
            haplo = haplo.split("|")[0]
            df = haps[haps["nomenclature"].isin(haplo.split("+"))]

            df = df.drop(
                columns=["haplotype", "nomenclature", "Jiang 2012 freq", "structure"]
            )
            df = df.set_axis(map(lambda i: f"KIR{i}", df.columns), axis=1)

            # cn
            name_id = self.getID(name)
            cn[name_id] = dict(df.sum(axis=0))

            # allele
            alleles = []
            for gene_name, gene_cn in cn[name_id].items():
                alleles.extend([gene_name] * int(gene_cn))
            guess_allele.append(
                {
                    "id": name_id,
                    "alleles": alleles,
                    "name": name,
                }
            )
        assert cn

        # cn
        df_cn = pd.DataFrame(cn)
        df_cn = df_cn.reset_index().rename(columns={"index": "gene"})
        logger.debug(f"[KPI] CN {df_cn}")
        df_cn.to_csv(output_name_cn + ".csv", index=False)

        # allele
        self.savePredictedAllele(guess_allele, output_name)
        return output_name

    def runAll(self, input_name: str) -> str:
        """Run all the script(Don't use this when building pipeline"""
        logger.info("[KPI] Download reference and code")
        index = self.download()
        samples = input_name
        logger.info(f"[KPI] Run nextflow {samples}")
        samples = self.run(samples, index=index)
        samples = self.mergeResult(samples, index=index)
        return samples
