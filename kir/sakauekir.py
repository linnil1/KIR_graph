"""The pipeline is modified from https://github.com/saorisakaue/KIR_project"""
import gzip
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .sakauekir_cn import getPloidy
from .kir_pipe import KirPipe


class SakaueKir(KirPipe):
    """Sakaue's KIR pipeline (rewritten by linnil1)"""

    name = "sakauekir"

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)
        self.images = {
            "bwa":         "quay.io/biocontainers/bwa:0.7.17--hed695b0_7",
            "gatk3":       "docker.io/broadinstitute/gatk3:3.6-0",
            "picard":      "quay.io/biocontainers/picard:2.27.3--hdfd78af_0",
            "samtools":    "quay.io/biocontainers/samtools:1.15.1--h1170115_0",
            "deepvariant": "docker.io/google/deepvariant:1.4.0",
        }
        self.version = "v1.0.0"

    def download(self, folder_base: str = "") -> str:
        """Download Database"""
        if folder_base and not folder_base.endswith("/"):
            folder_base += "/"
        folder = folder_base + "sakauekir_" + self.escapeName(self.version)
        if Path(folder).exists():
            return folder
        self.runShell(f"git clone https://github.com/saorisakaue/KIR_project {folder}")
        self.runShell(f"git checkout {self.version}", cwd=folder)
        return folder

    def bwa(self, input_name: str, index: str) -> str:
        """Step01: Align by bwa mem"""
        output_name = input_name + "." + self.escapeName(index) + ".bwa"
        if Path(output_name + ".bam").exists():
            return output_name
        f1, f2 = f"{input_name}.read.1.fq", f"{input_name}.read.2.fq"
        if not Path(f1).exists():
            f1, f2 = f"{input_name}.read.1.fq.gz", f"{input_name}.read.2.fq.gz"
        name_id = self.getID(input_name)
        rg = "@RG\\tID:" + name_id + "\\tSM:" + name_id
        self.runDocker(
            "bwa",
            f""" \
            bwa mem -t {self.getThreads()} \
                {index}/REF/KIR_seq_ref \
                -R "{rg}" \
                {f1} {f2} -o {output_name}.sam
            """,
        )
        self.samtobam(output_name)
        return output_name

    def addGroup(self, input_name: str) -> str:
        """Step01: Add readgroup"""
        output_name = input_name + ".rg"
        if Path(output_name + ".bam").exists():
            return output_name
        name_id = self.getID(input_name)
        self.runDocker(
            "picard",
            f""" \
            picard AddOrReplaceReadGroups \
                   I={input_name}.bam \
                   O={output_name}.bam \
                   RGLB={name_id} RGPL=ILLUMINA RGPU={name_id} RGSM={name_id} RGID={name_id} \
                   VALIDATION_STRINGENCY=LENIENT
            """,
        )
        return output_name

    def markDupliate(self, input_name: str) -> str:
        """Step01: MarkDuplicates"""
        output_name = input_name + ".md"
        if Path(output_name + ".bam").exists():
            return output_name
        self.runDocker(
            "picard",
            f""" \
            picard MarkDuplicates \
                   I={input_name}.bam \
                   O={output_name}.bam \
                   ASSUME_SORTED=false \
                   REMOVE_DUPLICATES=false \
                   CREATE_INDEX=True \
                   VALIDATION_STRINGENCY=LENIENT \
                   M={output_name}.metrics
            """,
        )
        return output_name

    def analysisTK(self, input_name: str, index: str) -> str:
        """Step01: Calculate some regions' depth"""
        output_name = input_name + ".coverage"
        if Path(output_name + ".vcf").exists():
            return output_name
        self.runDocker(
            "gatk3",
            f""" \
            java -jar /usr/GenomeAnalysisTK.jar \
                 -T DiagnoseTargets \
                 -I {input_name}.bam \
                 -o {output_name}.vcf \
                 -R {index}/REF/KIR_seq_ref.fasta \
                 -L {index}/REF/KIR_seq_ref.intervals
            """,
        )
        return output_name

    def getCoverage(self, input_name: str) -> str:
        """Step01.5: linnil1 Calculate the gene's depth before Step02"""
        output_name = input_name + ".depth"
        if Path(output_name + ".csv").exists():
            return output_name
        regions = []
        with open(input_name + ".vcf") as f_vcf:
            for line in f_vcf:
                if not line or line.startswith("#"):
                    continue
                # KIR2DL1	1	.	C	<DT>	.	PASS	END=159;IDP=17.63;IGC=0.604	IDP:LL:ZL	17.63:0:0
                fields = line.split("\t")
                info = dict(map(lambda i: i.split("="), fields[7].split(";")))
                gene = fields[0]
                length = float(info["END"]) - float(fields[1])
                regions.append(
                    {"gene": gene, "depth": float(info["IDP"]), "length": length}
                )

        # Paper: average depth with weighted by length
        df = pd.DataFrame(regions)
        df = df.groupby("gene").apply(lambda i: np.average(i.depth, weights=i.length))  # type: ignore
        df = df.reset_index()
        print(df)
        df.to_csv(output_name + ".csv", index=False, header=None)  # type: ignore
        return output_name

    def ploidyEstimate(self, input_name: str) -> str:
        """Step02: Determine copy number (with linnil1 modification python code)"""
        names = self.listFiles(input_name)
        output_base = self.replaceWildcard(input_name, "_merge_depth")
        output_name = output_base + ".ploidy"
        if Path(output_name + ".csv").exists():
            return output_name

        # merge df
        dfs = []
        for name in names:
            df = pd.read_csv(name + ".csv", header=None, index_col=0)
            df = df.set_axis([self.getID(name)], axis=1)
            dfs.append(df)
        df = pd.concat(dfs, axis=1)

        # Paper: normalize by 3DL3
        df = df / df.loc["KIR3DL3", :]  # type: ignore
        print(df)

        # Run my 02_determine_copy_number modified code
        ploidy = getPloidy(df, output_base)
        ploidy.loc["KIR3DL3", :] = 2
        ploidy = ploidy.fillna(0).astype(int)
        print(ploidy)

        # save
        ploidy.to_csv(output_name + ".csv")
        return output_name

    def beforeHC(self, input_name: str, ploidy_name: str) -> str:
        """Stpe02.5: Perpare copy number and bam data for Step03"""
        name_id = self.getID(input_name)
        output_name = (
            input_name
            + ".ploidy_"
            + self.escapeName(ploidy_name.format("same"))
            + ".gene.{}"
        )
        if Path(output_name.format("KIR3DL3") + ".json").exists():
            return output_name

        ploidy = pd.read_csv(ploidy_name.format(name_id) + ".csv", index_col=0)
        sample_ploidy = ploidy[name_id]
        print(sample_ploidy)

        for gene, ploidy in sample_ploidy.iteritems():
            if not ploidy:
                continue
            if gene == "KIR3DP1":
                continue
            with open(
                output_name.format(self.renameGene(str(gene))) + ".json", "w"
            ) as f_out:
                json.dump(
                    {
                        "id": name_id,
                        "gene": gene,
                        "input_name": input_name,
                        "bam": input_name + ".bam",
                        "ploidy": ploidy,
                    },
                    f_out,
                )
        return output_name

    def renameGene(self, gene: str) -> str:
        """
        Sometime they change the gene format.
        e.g.  "KIR2DL5A;KIR2DL5B" -> KIR2DL5AB
        """
        if "KIR2DL5A;KIR2DL5B" == gene:
            return "KIR2DL5AB"
        if "KIR2DS3;KIR2DS5" == gene:
            return "KIR2DS35"
        return gene

    def haplotypeCaller(self, input_name: str, index: str) -> str:
        """Step03: Haplotype caller"""
        output_name = input_name + ".hc"
        with open(input_name + ".json") as f:
            data = json.load(f)
        if Path(output_name + ".g.vcf.gz").exists():
            return output_name

        print(input_name, data)
        self.runDocker(
            "gatk3",
            f""" \
            java -jar /usr/GenomeAnalysisTK.jar \
                 -T HaplotypeCaller \
                 -I {data['bam']} \
                 -o {output_name}.g.vcf.gz \
                 -nct 2 \
                 -ploidy {data['ploidy']} \
                 -R {index}/REF/KIR_seq_ref.fasta \
                 -L '{index}/REF/{self.renameGene(data['gene'])}.intervals' \
                 --emitRefConfidence GVCF
            """,
        )
        return output_name

    def jointGenotype(self, input_name: str, index: str) -> str:
        """Step04: join all VCFs"""
        output_name = self.replaceWildcard(input_name, "_mergevcf") + ".gt"
        if Path(output_name + ".g.vcf.gz").exists():
            return output_name
        files = [name + ".g.vcf.gz" for name in self.listFiles(input_name)]
        self.runDocker(
            "gatk3",
            f""" \
            java -jar /usr/GenomeAnalysisTK.jar \
                 -T GenotypeGVCFs \
                 -R {index}/REF/KIR_seq_ref.fasta \
                 -allSites \
                 -o {output_name}.g.vcf.gz \
                 {" ".join(map(lambda i: "--variant " + i, files))}
        """,
        )
        return output_name

    def beforeCalling(self, input_name: str) -> str:
        """Step06.00: Call genotype"""
        output_name = input_name + ".genecall.{}"
        genes = [
            "KIR2DL1",
            "KIR2DL2",
            "KIR2DL3",
            "KIR2DL5A;KIR2DL5B",
            "KIR2DS1",
            "KIR2DS2",
            "KIR2DS3;KIR2DS5",
            "KIR2DS4",
            "KIR3DL1",
            "KIR3DL2",
            "KIR3DL3",
            "KIR3DS1",
            "KIR2DL5A;KIR2DL5B",
            "KIR2DS3;KIR2DS5",
            "KIR2DL4",
        ]
        if Path(output_name.format("KIR2DL4") + ".json").exists():
            return output_name
        for gene in sorted(set(genes)):
            with open(output_name.format(self.renameGene(gene)) + ".json", "w") as f:
                json.dump(
                    {
                        "input_name": input_name,
                        "vcf": input_name + ".g.vcf.gz",
                        "gene": gene,
                    },
                    f,
                )
        return output_name

    def calling(self, input_name: str, index: str) -> str:
        """Step06: Call genotype"""
        output_name = input_name
        if Path(output_name + ".alleles.tsv").exists():
            return output_name + ".alleles"

        # gene
        with open(input_name + ".json") as f:
            data = json.load(f)

        # get id from vcf
        with gzip.open(data["input_name"] + ".g.vcf.gz", "rt") as f:
            for i in f:
                i = str(i)
                if i.startswith("#CHROM"):
                    samples = [i.strip() for i in i.split("FORMAT")[-1].split()]
        assert samples
        name_id = samples[0]

        # call sakaue python by directed called
        self.runShell(
            f"""
            python kir/sakauekir_call.py \
                   {data['vcf']} \
                   '{index}/data/{data['gene']}.difpos.all.txt' \
                   {output_name}.dosage.tsv \
                   {output_name}.reference.tsv \
                   {output_name}.alleles.tsv \
                   '{data['gene']}' {name_id}
            """
        )
        return output_name + ".alleles"

    def mergeCalling(self, input_name: str) -> str:
        """Step06.2: Merge call results"""
        output_name = self.replaceWildcard(input_name, "_merge")
        if Path(output_name + ".tsv").exists():
            return output_name
        files = [i + ".tsv" for i in self.listFiles(input_name)]
        self.runShell(f"cat {' '.join(files)} > {output_name}.tsv")
        return output_name

    def readResult(
        self, filename: str, select_all: bool = False
    ) -> tuple[str, list[str]]:
        """linnil1 function: read alleles from alleles.tsv"""
        raw_df = pd.read_csv(filename, header=None, sep="\t", dtype=str)
        raw_df.columns = ["id", "gene", "alleles", "type"]  # type: ignore
        alleles = []
        name_id = ""
        for i in raw_df.itertuples():
            name_id = i.id
            if i.type == "known":
                alleles_text = i.alleles.replace("_", "*")
                possible_set = alleles_text.split("-or-")
            elif i.type == "potentially_novel":
                alleles_text = (
                    i.alleles.replace("Close_to_", "").replace("_", "*").split("[")[0]
                )
                possible_set = alleles_text.split("-OR-")
            else:
                raise ValueError(f"{i.type} type not found")

            if select_all:
                alleles.extend(
                    [
                        k
                        for i in possible_set
                        for j in i.split("/")
                        for k in j.split("-")
                    ]
                )
            else:
                alleles.extend([i.split("-")[0] for i in possible_set[0].split("/")])
        return name_id, alleles

    def mergeResult(self, input_name: str, select_all: bool = False) -> str:
        """Step06.3: Read result to our format"""
        if select_all:
            output_name = self.replaceWildcard(input_name, "_merge_called_full")
        else:
            output_name = self.replaceWildcard(input_name, "_merge_called")

        # if Path(output_name + ".tsv").exists():
        #     return output_name

        predict_list = []
        for name in self.listFiles(input_name):
            name_id, alleles = self.readResult(name + ".tsv", select_all=select_all)
            print(name_id, alleles)
            predict_list.append(
                {
                    "id": name_id,
                    "alleles": alleles,
                    "name": input_name.format(name_id),
                }
            )
        self.savePredictedAllele(predict_list, output_name)
        return output_name

    def deepVariant(self, input_name: str, index: str) -> str:
        """Step05: Deep variant for Structual variant calling.
        But I don't know how to merge 22bp record into main vcf"""
        output_name = input_name + ".dv"
        if Path(output_name + ".g.vcf.gz").exists():
            return output_name

        self.runDocker(
            "deepvariant",
            f""" \
            /opt/deepvariant/bin/run_deepvariant \
                --model_type=WGS \
                --ref {index}/REF/KIR_seq_ref.fasta \
                --reads {input_name}.bam \
                --output_vcf={output_name}.vcf.gz \
                --output_gvcf={output_name}.g.vcf.gz \
            """,
        )
        return output_name

    def runAll(self, input_name: str) -> str:
        """Run all the script(Don't use this when building pipeline"""
        folder = self.download()
        sample_bam = []
        samples = input_name
        for sample in self.listFiles   (samples):
            sample =  self.bwa         (sample, index=folder)
            sample =  self.addGroup    (sample)
            sample =  self.markDupliate(sample)
            sample_bam.append          (sample)
            sample =  self.analysisTK  (sample, index=folder)
            sample =  self.getCoverage (sample)

        samples += ".sakauekir_v1_0_0.bwa.rg.md"
        samples_cn = samples + ".coverage.depth"
        samples_cn = self.ploidyEstimate(samples_cn)
        for sample in sample_bam:
            samples_gene = self.beforeHC       (sample, samples_cn)
            sample       = self.deepVariant    (sample,               index=folder)
            for sample in  self.listFiles      (samples_gene):
                sample   = self.haplotypeCaller(sample,               index=folder)
            sample       = self.jointGenotype  (samples_gene + ".hc", index=folder)
            samples_call = self.beforeCalling  (sample)
            for sample in  self.listFiles      (samples_call):
                sample   = self.calling        (sample,               index=folder)
            sample       = self.mergeCalling   (samples_call + ".alleles")

        samples = (
            samples
            + ".ploidy_example_data_test_merge_depth_sakauekir_v1_0_0_bwa_rg_md_coverage_depth_ploidy"
            + ".gene_mergevcf.hc.gt.genecall_merge.alleles"
        )
        sample = self.mergeResult(samples)
        return sample
