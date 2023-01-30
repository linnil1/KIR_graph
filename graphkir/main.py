"""
GRAPH-KIR: Main function

TODO: update after kg_main OK
"""
from glob import glob
from pathlib import Path
import argparse

import pandas as pd
import plotly.graph_objects as go

# order of pipeline
from .kir_msa import buildKirMsa
from .msa_leftalign import genemsaLeftAlign
from .msa2hisat import msa2HisatReference, buildHisatIndex
from .wgs import downloadHg19, bwa, bwaIndex, extractFromHg19, bam2fastq
from .hisat2 import hisatMap, extractVariantFromBam, readExons
from .kir_cn import predictSamplesCN, loadCN, filterDepth, bam2Depth
from .kir_typing import selectKirTypingModel
from .utils import setThreads, getThreads, mergeCN, mergeAllele, setEngine
from .plot import plotCN, plotReadMappingStat, showPlot, plotGeneDepths


def buildMSA(
    msa_type: str = "ab_2dl1s1",
    index_folder: str = "index",
    threads: int = 1,
    add_exon_only_sequences: bool = True,
    ipd_version: str = "2100",
) -> str:
    Path(index_folder).mkdir(exist_ok=True)
    """Build MSA from source"""
    if not add_exon_only_sequences:
        msa_index = f"{index_folder}/kir_{ipd_version}_{msa_type}"
        if not glob(msa_index + ".KIR*"):
            buildKirMsa(
                msa_type, msa_index,
                version=ipd_version,
                threads=threads,
            )
    else:
        exon_name = f"{index_folder}/kir_{ipd_version}_withexon"
        msa_index = exon_name + "_" + msa_type
        if not glob(exon_name + "_ab.KIR*"):
            buildKirMsa(
                "ab", exon_name + "_ab",
                version=ipd_version,
                full_length_only=False,
                threads=threads,
            )
        if not glob(msa_index + ".KIR*"):
            buildKirMsa(
                msa_type, msa_index,
                version=ipd_version,
                input_msa_prefix=exon_name + "_ab",
                threads=threads,
            )

    msa_index1 = msa_index + ".leftalign"
    if not glob(msa_index1 + ".KIR*"):
        genemsaLeftAlign(msa_index, msa_index1)
    msa_index = msa_index1
    return msa_index


def buildIndex(
    msa_type: str = "ab_2dl1s1",
    index_folder: str = "index",
    threads: int = 1,
    add_exon_only_sequences: bool = True,
    ipd_version: str = "2100",
) -> tuple[str, str]:
    """
    Build graph data

    * MSA
    * HISAT2 formatted files (for hisat-build)
    * HISAT2 index

    Returns:
      HISAT2 files prefix
      HISAT2 index prefix
    """
    if not add_exon_only_sequences:
        msa_index = f"{index_folder}/kir_2100_{msa_type}.leftalign"
    else:
        msa_index = f"{index_folder}/kir_2100_withexon_{msa_type}.leftalign"

    ref_index = msa_index + ".mut01"
    if not Path(ref_index + ".haplotype").exists():
        if not glob(msa_index + ".KIR*"):
            buildMSA(
                msa_type, index_folder,
                threads=threads,
                add_exon_only_sequences=add_exon_only_sequences,
                ipd_version=ipd_version,
            )
        msa2HisatReference(msa_index, ref_index)

    index = ref_index + ".graph"
    if not Path(index + ".8.ht2").exists():
        buildHisatIndex(ref_index, index)

    return ref_index, index


def buildGenomeIndex(index_folder: str = "index") -> str:
    """Download hs37d5.fa and build hs37d5 bwa index"""
    Path(index_folder).mkdir(exist_ok=True)
    wgs_index = f"{index_folder}/hs37d5.fa.gz"
    if not Path(wgs_index).exists():
        wgs_index = downloadHg19(index_folder)
    if not Path(wgs_index + ".bwt").exists():
        bwaIndex(wgs_index, wgs_index)
    return wgs_index


def getCommonName(r1: str, r2: str) -> str:
    """
    Extract the name for r1 and r2.

    The name is the longest common string between r1 and r2.

    Example:
      Input:
      * "data/linnil1_syn_30x.00.read.r1.fq" "data/linnil1_syn_30x.00.read.r2.fq"
      * "data/linnil1_syn_30x.01.read.r1.fq" "data/linnil1_syn_30x.01.read.r2.fq"

      Output's name:
      * "data/linnil1_syn_30x.00"
      * "data/linnil1_syn_30x.01"

    Returns:
      list of tuple[name, r1, r2]
    """
    name = ""
    for s1, s2 in zip(r1.split("."), r2.split(".")):
        if s1 == s2:
            if name:
                name += "." + s1
            else:
                name = s1
        else:
            return name
    return name


def replaceParentFolder(filename: str, new_folder: str) -> str:
    """a/b/c.tsv -> new_folder/c.tsv"""
    return str(Path(new_folder) / Path(filename).name)


def createParser() -> argparse.ArgumentParser:
    """Setup arguments"""
    parser = argparse.ArgumentParser(
        description="Run Graph KIR",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--thread",
        default=1,
        help="Number of threads",
    )
    parser.add_argument(
        "--engine",
        default="podman",
        help="podman,docker,singularity,local(samtools,hisat2... installed)",
    )
    parser.add_argument(
        "--ipd-version",
        default="2100",
        help="IPD*KIR version. Branch name in https://github.com/ANHIG/IPDKIR are all available.",
    )
    parser.add_argument(
        "--r1",
        action="append",
        help="Read 1 fastq",
    )
    parser.add_argument(
        "--r2",
        action="append",
        help="Read 2 fastq",
    )
    parser.add_argument(
        "--input-csv",
        help="All data are written in csv. Including id, fastq, cnfile",
    )
    parser.add_argument(
        "--index-folder",
        default="index",
        help="The path to index",
    )
    parser.add_argument(
        "--output-folder",
        help="The path for saving data. Otherwise save in the same folder with reads.",
    )
    parser.add_argument(
        "--output-cohort-name",
        help="Cohort name (Default: {folder}/cohort.xx)",
    )
    parser.add_argument(
        "--msa-type",
        default="ab_2dl1s1",
        help="Type of MSA: merge, split, ab, ab_2dl1s1",
    )
    parser.add_argument(
        "--msa-no-exon-only-allele",
        action="store_true",
        help="Do not add exon-only alleles into allele set",
    )
    parser.add_argument(
        "--cn-exon",
        action="store_true",
        help="Select exon-only region of CN prediction",
    )
    parser.add_argument(
        "--cn-multi",
        action="store_true",
        help="Do not remove multiple aligned reads for predicting CN",
    )
    parser.add_argument(
        "--cn-individually",
        action="store_true",
        help="Predict CN separtely instead of cohort-wisely",
    )
    parser.add_argument(
        "--cn-select",
        default="p75",
        help="Select 75 percentile(p75), mean(mean), median(median) of depths",
    )
    parser.add_argument(
        "--cn-cluster",
        default="CNgroup",
        help="CN prediction model (CNgroup, kde)",
    )
    parser.add_argument(
        "--cn-provided",
        nargs="*",
        help="Provided CN prediction in TSV",
    )
    parser.add_argument(
        "--allele-method",
        default="pv",
        help="Max Likelihood by variants(pv), report(report)",
    )
    parser.add_argument(
        "--allele-no-exon-only-allele",
        action="store_true",
        help="Do not use exon-first algorithm",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Plot some intermediate result.",
    )
    parser.add_argument(
        "--extract-skip",
        action="store_true",
        help="Filtering KIR reads from the fastq",
    )
    parser.add_argument(
        "--extract-index",
        default="",
        help="A bwa-indexed hs19 index path (Default=index_folder/hs37d5.fa.gz)",
    )
    return parser


def main(args: argparse.Namespace) -> list[go.Figure]:
    """
    * Read Mapping
    * Extract variants
    * CN prediction
    * allele typing
    * merge
    """
    # Configuration
    setThreads(args.thread)
    setEngine(args.engine)

    # Init
    cn_files = []
    if not args.input_csv:
        assert len(args.r1)
        assert len(args.r1) == len(args.r2)
        reads = list(zip(args.r1, args.r2))
        names = list(map(lambda i: getCommonName(i[0], i[1]), reads))
        if args.cn_provided:
            cn_files = args.cn_provided
        else:
            cn_files = [""] * len(names)
    else:
        df = pd.read_csv(args.input_csv)
        names = list(df["name"])
        reads = list(zip(df["r1"], df["r2"]))
        if "cnfile" in df.columns:
            # Doesn't require all cnfile existed
            cn_files = list(df["cnfile"].fillna(""))
    assert len(cn_files) == len(names)

    if args.output_folder:
        output_folder = args.output_folder
        Path(output_folder).mkdir(exist_ok=True)
        names = list(map(lambda i: replaceParentFolder(i, output_folder), names))
    else:
        output_folder = str(Path(names[0]).parent)

    if args.output_cohort_name:
        cohort_name = args.output_cohort_name
        Path(cohort_name).parent.mkdir(exist_ok=True)
    else:
        cohort_name = str(Path(output_folder) / "cohort")
    bam_files = []
    processed_reads = []
    depth_files = []
    allele_files = []
    figs = []

    # extraction step (optional)
    if not args.extract_skip:
        # Prepare wgs Index
        if not args.extract_index:
            wgs_index = buildGenomeIndex(args.index_folder)
        else:
            wgs_index = args.extract_index

        # read mapping and extract to fastq
        new_names = []
        new_reads = []
        for name, (fq1, fq2) in zip(names, reads):
            suffix = "." + wgs_index.replace(".", "_").replace("/", "_")
            bwa(wgs_index, fq1, fq2, name + suffix, threads=getThreads())
            name += suffix

            # extract
            suffix = ".extract"
            extractFromHg19(
                name + ".bam", name + suffix, "hs37d5", threads=getThreads()
            )
            name += suffix
            bam2fastq(name + ".bam", name, threads=getThreads())
            new_names.append(name)
            new_reads.append((name + ".read.1.fq.gz", name + ".read.2.fq.gz"))
        reads = new_reads
        names = new_names

    # Prepare Index
    ref_index, index = buildIndex(
        args.msa_type, args.index_folder,
        threads=getThreads(),
        add_exon_only_sequences=not args.msa_no_exon_only_allele,
        ipd_version=args.ipd_version,
    )

    # Read Mapping and filtering
    for name, (fq1, fq2) in zip(names, reads):
        suffix = "." + index.replace(".", "_").replace("/", "_")
        hisatMap(index, fq1, fq2, name + suffix + ".bam", threads=getThreads())
        name += suffix
        bam_files.append(name + ".bam")

        suffix = ".variant"
        extractVariantFromBam(
            ref_index, name + ".bam", name + suffix, error_correction=False
        )
        name += suffix
        processed_reads.append(name)  # it has multiple format

        # read depth pre-process
        name += ".no_multi" if not args.cn_multi else ""
        suffix = ".depth"
        bam2Depth(name + ".bam", name + suffix + ".tsv")
        name += suffix

        # extract exon regions
        if args.cn_exon:
            exon_regions = readExons(ref_index)
            suffix = ".exon"
            filterDepth(name + ".tsv", name + suffix + ".tsv", exon_regions)
            name += suffix
        depth_files.append(name + ".tsv")
    figs.extend(plotReadMappingStat(bam_files, [fq1 for fq1, _ in reads]))

    # Copy Number determination
    if all(cn_files):
        pass
    elif not args.cn_individually:
        for i, depth_file in enumerate(depth_files):
            if cn_files[i]:  # CN file exists, skip
                continue
            suffix = f".{args.cn_select}.{args.cn_cluster}"
            name = str(Path(depth_file).with_suffix(suffix))
            predictSamplesCN(
                [depth_file],
                [name + ".tsv"],
                cluster_method=args.cn_cluster,
                assume_3DL3_diploid=True,
                save_cn_model_path=name + ".json",
                select_mode=args.cn_select,
            )
            cn_files[i] = name + ".tsv"
            figs.extend(plotGeneDepths(depth_file))
            figs.extend(plotCN(name + ".json"))
    else:
        suffix = f".{args.cn_select}.cohort.{args.cn_cluster}"
        cn_cohort_name = cohort_name + suffix
        cn_files = [str(Path(path).with_suffix(suffix + ".tsv")) for path in depth_files]
        predictSamplesCN(
            [depth_files[i] for i, cnf in enumerate(cn_files) if cnf],
            [cnf            for i, cnf in enumerate(cn_files) if cnf],
            save_cn_model_path=cn_cohort_name + ".json",
            cluster_method=args.cn_cluster,
            select_mode=args.cn_select,
        )
        figs.append(plotCN(cn_cohort_name + ".json"))
    print(f"Saved in {cohort_name}.cn.tsv")
    print(mergeCN(cn_files, cohort_name + ".cn.tsv"))

    # Allele Typing
    for name, cn_file in zip(processed_reads, cn_files):
        # suffix = ".cn_" + cn_file.replace("/", "_").replace(".", "_") + "."
        # filename too long
        suffix = (
            ".cn"
            + cn_file[len(getCommonName(name, cn_file)) :]
            .replace("/", "_")
            .replace(".", "_")
            + "."
        )
        method = args.allele_method
        if not args.allele_no_exon_only_allele:
            assert method == "pv"
            method += "_exonfirst_1"  #_1.2
        suffix += method
        t = selectKirTypingModel(
            method, name + ".json",
            top_n=600,
            variant_correction=True,
        )
        cn = loadCN(cn_file)
        called_alleles = t.typing(cn)
        name += suffix
        # t.save(name + ".json")
        df = pd.DataFrame(
            {
                "name": [name],
                "alleles": ["_".join(called_alleles)],
            }
        )
        df.to_csv(name + ".tsv", sep="\t", index=False)
        allele_files.append(name + ".tsv")
    print(f"Saved in {cohort_name}.allele.tsv")
    print(mergeAllele(allele_files, cohort_name + ".allele.tsv"))
    return figs


def entrypoint() -> None:
    """Command line entry function"""
    parser = createParser()
    args = parser.parse_args()
    figs = main(args)
    if args.plot:
        showPlot(figs)


if __name__ == "__main__":
    entrypoint()
