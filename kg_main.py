import os
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

import kg_build_index
import kg_typping

from kg_build_msa import (
    kirToMultiMsa,
    kirToSingleMsa,
    kirMerge2dl1s1,
)

from kg_utils import (
    runShell,
    runDocker,
    runAll,
    getSamples,
    threads,
    samtobam,
)


def linkSamples():
    os.makedirs("data", exist_ok=True)
    index = "data/linnil1_syn_wide"
    if getSamples(index, ".fq", strict=False):
        return index

    files = getSamples("linnil1_syn_wide/linnil1_syn_wide", ".fq", strict=False, return_name=True)
    for f, id in files:
        runShell(f"ln -s ../{f} {index}.{id}.fq")
    return index


def link10Samples(sample_index):
    suffix = ".test10"
    if getSamples(sample_index + suffix, ".fq", strict=False):
        return sample_index + suffix

    # yes, brute force
    for id in range(10):
        runShell(f"ln -s {os.path.basename(sample_index)}.{id:02d}.read.1.fq {sample_index}{suffix}.{id:02d}.read.1.fq")
        runShell(f"ln -s {os.path.basename(sample_index)}.{id:02d}.read.2.fq {sample_index}{suffix}.{id:02d}.read.2.fq")
    return sample_index + suffix


def link30xSamples():
    os.makedirs("data", exist_ok=True)
    index = "data/linnil1_syn_30x"
    if getSamples(index, ".fq", strict=False):
        return index

    files = getSamples("linnil1_syn_30x/linnil1_syn_30x", ".fq", strict=False, return_name=True)
    for f, id in files:
        runShell(f"ln -s ../{f} {index}.{id}.fq")
    return index


def hisatMapAll(index, sample_index, suffix=""):
    new_suffix = ""
    for f in getSamples(sample_index, suffix + ".read.1.fq"):
        name = f[:-len(".read.1.fq")]
        new_suffix = hisatMap(index, name)
    return new_suffix


def hisatMap(index, name):
    suffix = "." + index.replace("/", "_")
    if os.path.exists(f"{name}{suffix}.bam"):
        return suffix
    f1, f2 = name + ".read.1.fq", name + ".read.2.fq"
    runDocker("hisat", f"""\
              hisat2 --threads {threads} -x {index}.graph -1 {f1} -2 {f2} \
              --no-spliced-alignment --max-altstried 64 --haplotype \
              -S {name}{suffix}.sam """)
    samtobam(f"{name}{suffix}")
    return suffix


def bowtie2BuildFull(index, from_index="index/kir_2100_raw.mut01"):
    # No matter what index, this file will have the same sequences
    index = f"{index}/kir_2100_raw_full"
    if getSamples(index, strict=False):
        return index
    old_index = from_index
    runDocker("bowtie",
              f"bowtie2-build {old_index}_sequences.fa {index} --threads {threads}")
    return index


def bowtie2BuildConsensus(index, from_index="index/kir_2100_raw.mut01"):
    index = f"{index}/kir_2100_raw_cons"
    if getSamples(index, strict=False):
        return index
    old_index = from_index
    runDocker("bowtie",
              f"bowtie2-build {old_index}_backbone.fa {index} --threads {threads}")
    return index


def bowtie2All(index, sample_index, suffix="", use_arg="default"):
    def runBowtie(name):
        name = name[:-len(".read.1.fq")]
        new_suffix = bowtie2(index, name, use_arg)
        return new_suffix

    new_suffix = ""
    for i in runAll(runBowtie,
                    getSamples(sample_index, suffix + ".read.1.fq"),
                    concurrent=False):
        new_suffix = i
    return new_suffix


def bowtie2(index, name, use_arg="default"):
    suffix = "." + index.replace("/", "_") + ".bowtie2"
    args = ""
    if use_arg == "ping":
        args = "  -5 0  -3 6  -N 0  --end-to-end  --score-min L,-2,-0.08  " + \
               "  -I 75  -X 1000  -a  --np 1  --mp 2,2  --rdg 1,1  --rfg 1,1  "
        suffix += "_ping"
    if os.path.exists(f"{name}{suffix}.bam"):
        return suffix

    # main
    f1, f2 = name + ".read.1.fq", name + ".read.2.fq"
    runDocker("bowtie",
              f"bowtie2 {args} --threads {threads} -x {index} "
              f"-1 {f1} -2 {f2} -a -S {name}{suffix}.sam")
    samtobam(name + suffix)
    return suffix


def hisatTyping(index, sample_index, suffix=""):
    new_suffix = ""
    kg_typping.setIndex(index)
    samples = getSamples(sample_index, suffix + ".bam")
    for i in runAll(kg_typping.main, samples):
        new_suffix = i
    return new_suffix


if __name__ == "__main__":
    # sample_index = linkSamples()
    # sample_index = link10Samples(sample_index)
    sample_index = link30xSamples()

    index = "index"
    suffix = ""
    os.makedirs("index", exist_ok=True)
    # Using IPDKIR MSA
    # index = kirToMultiMsa(index)           # index = "index/kir_2100_raw.mut01"
    # Merge 2DL1 2DS1
    index = kirMerge2dl1s1(index)          # index = "index/kir_2100_2dl1s1.mut01"
    # Split 2DL5A and 2DL5B
    # index = kirToMultiMsa(split_2DL5=True) # index = "index/kir_2100_ab"
    # Merge all KIR
    # index = kirToSingleMsa(index)          # index = "index/kir_2100_merge.mut01"


    index = kg_build_index.main(index)
    suffix += hisatMapAll(index, sample_index, suffix)
    suffix += hisatTyping(index, sample_index, suffix)
    print(index, sample_index, suffix)


    # Different method mapping test
    index = "index"
    suffix = ""
    # index = bowtie2BuildConsensus(index, from_index="index/kir_2100_raw.mut01")  # index = "index/kir_2100_raw_cons"
    # index = bowtie2BuildFull(index, from_index="index/kir_2100_raw.mut01")       # index = "index/kir_2100_raw_full"

    # suffix += bowtie2All(index, sample_index, suffix)
    # suffix += bowtie2All(index, sample_index, suffix, use_arg="ping")
    """
    index = "index/kir_2100_muscle"
    # kirToSingleMsa(method='muscle')
    # bamFilter("-f 0x2", ".nosingle")
    # suffix += ".nosec"
    # bamFilter("-f 0x2 -F 256", ".sec")
    """
    runShell("stty echo opost")
