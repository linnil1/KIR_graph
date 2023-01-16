"""
Try extract KIR3DX1 -> to short reads -> mapped on graph -> stat

# UCSC hg19
# KIR3DX1                  chr19  55043909   55057024

# UCSC hg38
602   KIR3DX1                 chr19  54532692   54545741   13049
576   KIR3DX1  chr19_GL949746v1_alt    438900     451940   13040
569   KIR3DX1  chr19_GL949752v1_alt    437417     450458   13041
1051  KIR3DX1  chr19_KI270938v1_alt    515023     528138   13115
"""
from pathlib import Path
from functools import partial

from namepipe import compose

from kg_main import createSampleFastq, hisatMapWrap, extractVariant
from kg_utils import runDocker, runShell, back


def extractKIR3DX1(input_name, folder="", version="hg19"):
    """Extract KIR3DX1 and flanking = 500bp"""
    if not folder:
        output_name = input_name + ".3dx1_500bp"
    else:
        output_name = folder + "/" + Path(input_name).name + ".3dx1_f500bp"
    if Path(output_name + ".fa").exists():
        return output_name + ".{}"
    if version == "hg19":
        runShell(f"gunzip -k {input_name}.fa.gz")
        runDocker(
            "samtools",
            f"samtools faidx {input_name}.fa 19:55043409-55057524 -o {output_name}.19.fa",
        )
    elif version == "hg38":
        # runShell(f"gunzip -k {input_name}.fa.gz")
        runDocker(
            "samtools",
            f"samtools faidx {input_name}.fa"
            "   chr19:54532692-54545741"
            "   chr19_GL949746v1_alt:438400-452340"
            "   chr19_GL949752v1_alt:436917-450958"
            "   chr19_KI270938v1_alt:514523-528638"
            f"  -o {output_name}.38.fa",
        )
    else:
        raise ValueError("version not found")
    return output_name + ".{}"


def downloadUCSCHg38(folder):
    """Download UCSC hg38.p13 genome"""
    output_name = folder + "/hg38.ucsc"
    if Path(output_name + ".fa.gz").exists():
        return output_name
    runShell(
        f"wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz -O {output_name}.fa.gz"
    )
    return output_name


def showStat(input_name):
    """Display the samtoolss flagstat only"""
    runDocker("samtools", f"samtools flagstat {input_name}.bam")
    return input_name


if __name__ == "__main__":
    ref_index = "index/kir_2100_withexon_ab_2dl1s1.leftalign.mut01"
    index = ref_index + ".graph"
    version = "hg38"
    if version == "hg38":
        genome = downloadUCSCHg38("index")
    else:
        genome = ("index/hs37d5",)  # this file is already downloaded at main.py
    compose(
        [
            genome,
            partial(extractKIR3DX1, folder="data", version=version),
            partial(createSampleFastq, depth=30, seed=444),
            back,
            partial(hisatMapWrap, index=str(index)),
            showStat,
        ]
    )
