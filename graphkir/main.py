"""
GRAPH-KIR: Main function

TODO: update after kg_main OK
"""
from glob import glob
from pathlib import Path
import logging
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
from .utils import setThreads, getThreads, mergeCN, mergeAllele, setEngine, logger
from .plot import plotCN, plotReadMappingStat, showPlot, plotGeneDepths


def buildMSA(
    msa_type: str = "ab_2dl1s1",
    index_folder: str = "index",
    add_exon_only_sequences: bool = True,
    ipd_version: str = "2100",
) -> str:
    """Build MSA from source"""
    Path(index_folder).mkdir(exist_ok=True)
    if not add_exon_only_sequences:
        msa_index = f"{index_folder}/kir_{ipd_version}_{msa_type}"
        if not glob(msa_index + ".KIR*"):
            buildKirMsa(
                msa_type,
                msa_index,
                version=ipd_version,
                threads=getThreads(),
            )
    else:
        exon_name = f"{index_folder}/kir_{ipd_version}_withexon"
        msa_index = exon_name + "_" + msa_type
        if not glob(exon_name + "_ab.KIR*"):
            buildKirMsa(
                "ab",
                exon_name + "_ab",
                version=ipd_version,
                full_length_only=False,
                threads=getThreads(),
            )
        if not glob(msa_index + ".KIR*"):
            buildKirMsa(
                msa_type,
                msa_index,
                version=ipd_version,
                input_msa_prefix=exon_name + "_ab",
                threads=getThreads(),
            )

    msa_index1 = msa_index + ".leftalign"
    if not glob(msa_index1 + ".KIR*"):
        genemsaLeftAlign(msa_index, msa_index1)
    msa_index = msa_index1
    return msa_index


def buildGenomeIndex(index_folder: str = "index") -> str:
    """Download hs37d5.fa and build hs37d5 bwa index"""
    Path(index_folder).mkdir(exist_ok=True)
    wgs_index = f"{index_folder}/hs37d5.fa.gz"
    if not Path(wgs_index).exists():
        # logger.info(f"[WGS] Download hs37d5")
        wgs_index = downloadHg19(index_folder)
    if not Path(wgs_index + ".bwt").exists():
        logger.info(f"[WGS] Build {wgs_index} bwa index")
        bwaIndex(wgs_index, wgs_index)
    return wgs_index


def runWGS(
    names: list[str], reads: list[tuple[str, str]], index_wgs: str
) -> tuple[list[str], list[tuple[str, str]]]:
    new_names = []
    new_reads = []
    for name, (fq1, fq2) in zip(names, reads):
        # read mapping
        logger.info(f"[WGS] Run BWA on index {index_wgs} ({name})")
        suffix = "." + index_wgs.replace(".", "_").replace("/", "_")
        bwa(index_wgs, fq1, fq2, name + suffix, threads=getThreads())
        name += suffix

        # extract
        suffix = ".extract"
        # logger.info(f"[WGS] Extract chr19... from {input_bam}")  # written in func
        extractFromHg19(name + ".bam", name + suffix, "hs37d5", threads=getThreads())
        name += suffix
        logger.info(f"[WGS] Extract read from {name}.bam")
        bam2fastq(name + ".bam", name, threads=getThreads())
        new_names.append(name)
        new_reads.append((name + ".read.1.fq.gz", name + ".read.2.fq.gz"))
    return new_names, new_reads


def readMapping(
    names: list[str],
    reads: list[tuple[str, str]],
    index: str,
    index_ref: str,
    exon_region_only: bool = False,
) -> tuple[list[str], list[str], list[str]]:
    """
    1. Graph read mapping
    2. Read filtering(edit-distance thresholding no-mutliple-mapping)
    3. Calculate read depth
    """
    bam_files = []
    depth_files = []
    processed_bam = []
    for name, (fq1, fq2) in zip(names, reads):
        # mapping and filter
        suffix = "." + index.replace(".", "_").replace("/", "_")
        logger.info(f"[Graph] Run graph mapping on index {index} ({name})")
        hisatMap(index, fq1, fq2, name + suffix + ".bam", threads=getThreads())
        name += suffix
        bam_files.append(name + ".bam")

        suffix = ".variant"
        logger.info(f"[Graph] Filter mapping ({name})")
        extractVariantFromBam(
            index_ref, name + ".bam", name + suffix, error_correction=False
        )
        name += suffix
        processed_bam.append(name)  # it has multiple format

        # read depth pre-process
        name += ".no_multi"
        suffix = ".depth"
        logger.info(f"[Graph] Calculate read depth to {name}{suffix}.tsv")
        bam2Depth(name + ".bam", name + suffix + ".tsv")
        name += suffix

        # extract exon regions
        if exon_region_only:
            exon_regions = readExons(index_ref)
            suffix = ".exon"
            logger.info(f"[Graph] Filter exon read to {name}{suffix}.tsv")
            filterDepth(name + ".tsv", name + suffix + ".tsv", exon_regions)
            name += suffix
        depth_files.append(name + ".tsv")
    return bam_files, processed_bam, depth_files


def alleleTyping(
    processed_bam: list[str],
    cn_files: list[str],
    method: str = "full",
) -> list[str]:
    """Allele typing"""
    allele_files = []
    for name, cn_file in zip(processed_bam, cn_files):
        logger.debug(f"[Allele] Allele typing ({method}) with CN {cn_file} ({name})")
        # suffix = ".cn_" + cn_file.replace("/", "_").replace(".", "_") + "."
        # filename too long
        suffix = (
            ".cn"
            + cn_file[len(getCommonName(name, cn_file)) :]
            .replace("/", "_")
            .replace(".", "_")
            + "."
        )
        if method == "exonfirst":
            method += "_1"
        suffix += method
        t = selectKirTypingModel(
            method,
            name + ".json",
            top_n=600,
            variant_correction=True,
        )
        cn = loadCN(cn_file)
        called_alleles = t.typing(cn)
        logger.info(f"[Allele] {called_alleles} ({name})")
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

        logger.info(f"[Allele] All possible allele set in Allele typing saved in [name].possible.tsv")
        possible_list = t.getAllPossibleTyping()
        df_possible = pd.DataFrame(possible_list)
        df_possible = df_possible.fillna("")
        df_possible.to_csv(name + ".possible.tsv", index=False, sep="\t")
    return allele_files


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
        default="local",
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
        help="The path to graph index",
    )
    parser.add_argument(
        "--index-wgs",
        default="",
        help="A bwa-indexed hs37d5 index path (Default=index_folder/hs37d5.fa.gz)",
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
        choices=["merge", "split", "ab", "ab_2dl1s1"],
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
        "--cn-group-dev",
        default=0.08,
        help="Deviation of CNgroup""",
    )
    parser.add_argument(
        "--cn-group-3dl3-not-diploid",
        action="store_true",
        help="Not assume KIR3DL3 as diploid gene in CNgroup""",
    )
    parser.add_argument(
        "--cn-provided",
        nargs="*",
        help="Provided CN prediction in TSV",
    )
    parser.add_argument(
        "--allele-method",
        default="full",
        choices=["full", "exonfirst", "report"],
        help="Max Likelihood by full variants(full), "
        "exon variant first than full (exon-first) or hisat-report(report)",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Plot some intermediate result.",
    )
    parser.add_argument(
        "--step-skip-extraction",
        action="store_true",
        help="Skip filtering KIR reads from the (WGS) fastq",
    )
    parser.add_argument(
        "--step-skip-typing",
        action="store_true",
        help="Skip allele-typing",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=logging._nameToLevel.keys(),
        help="Set log level",
    )
    return parser


def main(args: argparse.Namespace) -> list[go.Figure]:
    """
    * WGS read mapping and KIR extraction
    * Graph Read Mapping
    * Extract variants
    * Calculate read depth + CN prediction
    * Allele typing
    * Merge result
    """
    # Configuration
    setThreads(args.thread)
    setEngine(args.engine)
    logging.basicConfig(level=args.log_level)
    logger.debug(f"[Main] {args}")

    # Reorganize input
    # name of sample, reads, and cn_file
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
        else:
            cn_files = [""] * len(names)
    if not names:
        logger.error(f"[Main] 0 Samples")
        exit()
    assert len(cn_files) == len(names)
    logger.info(f"[Main] Samples: {names}")

    # output name
    if args.output_folder:
        output_folder = args.output_folder
        Path(output_folder).mkdir(exist_ok=True)
        names = list(map(lambda i: replaceParentFolder(i, output_folder), names))
    else:
        output_folder = str(Path(names[0]).parent)

    # cohort name (final output name)
    if args.output_cohort_name:
        cohort_name = args.output_cohort_name
        Path(cohort_name).parent.mkdir(exist_ok=True)
    else:
        cohort_name = str(Path(output_folder) / "cohort")

    # init other
    allele_files = []
    figs = []

    # extraction step (optional)
    if not args.step_skip_extraction:
        # Prepare wgs Index
        if not args.index_wgs:
            index_wgs = buildGenomeIndex(args.index_folder)
        else:
            index_wgs = args.index_wgs

        # read mapping and extract to fastq
        names, reads = runWGS(names, reads, index_wgs)

    # Prepare MSA
    if args.msa_no_exon_only_allele:
        index_msa = (
            f"{args.index_folder}/kir_{args.ipd_version}_{args.msa_type}.leftalign"
        )
    else:
        index_msa = f"{args.index_folder}/kir_{args.ipd_version}_withexon_{args.msa_type}.leftalign"
    if not glob(index_msa + ".KIR*"):
        buildMSA(
            args.msa_type,
            args.index_folder,
            add_exon_only_sequences=not args.msa_no_exon_only_allele,
            ipd_version=args.ipd_version,
        )

    # MSA -> hisat
    index_ref = index_msa + ".mut01"
    if not Path(index_ref + ".haplotype").exists():
        logger.info(f"[MSA] MSA -> HISAT2 format")
        msa2HisatReference(index_msa, index_ref)

    # hisat -> graph index
    index = index_ref + ".graph"
    if not Path(index + ".8.ht2").exists():
        logger.info(f"[MSA] HISAT2 format -> HISAT2 graph index")
        buildHisatIndex(index_ref, index)

    # Read Mapping and filtering
    bam_files, processed_bam, depth_files = readMapping(
        names,
        reads,
        index,
        index_ref,
        exon_region_only=args.cn_exon,
    )
    figs.extend(plotReadMappingStat(bam_files, [fq1 for fq1, _ in reads]))

    cluster_method_kwargs = {
        "base_dev": float(args.cn_group_dev),
        "start_base": 2,
    }
    # Copy Number determination
    if all(cn_files):
        pass
    elif args.cn_individually:
        # per sample
        for i, depth_file in enumerate(depth_files):
            if cn_files[i]:  # CN file exists, skip
                continue
            suffix = f".{args.cn_select}.{args.cn_cluster}"
            name = str(Path(depth_file).with_suffix(suffix))
            logger.info(f"[CN] Copy number estimation per sample ({name})")
            predictSamplesCN(
                [depth_file],
                [name + ".tsv"],
                cluster_method=args.cn_cluster,
                cluster_method_kwargs=cluster_method_kwargs,
                assume_3DL3_diploid=not args.cn_group_3dl3_not_diploid,
                save_cn_model_path=name + ".json",
                select_mode=args.cn_select,
            )
            cn_files[i] = name + ".tsv"
            figs.extend(plotGeneDepths(depth_file))
            figs.extend(plotCN(name + ".json"))
    else:
        # by cohort
        suffix = f".{args.cn_select}.cohort.{args.cn_cluster}"
        cn_cohort_name = cohort_name + suffix
        cn_files = [str(Path(path).with_suffix(suffix + ".tsv")) for path in depth_files]
        logger.info(f"[CN] Copy number estimation by cohort ({cn_cohort_name})")
        predictSamplesCN(
            [depth_files[i] for i, cnf in enumerate(cn_files) if cnf],
            [cnf            for i, cnf in enumerate(cn_files) if cnf],
            cluster_method=args.cn_cluster,
            cluster_method_kwargs=cluster_method_kwargs,
            save_cn_model_path=cn_cohort_name + ".json",
            select_mode=args.cn_select,
        )
        figs.append(plotCN(cn_cohort_name + ".json"))
    logger.debug(f"[CN] Copy number files: {cn_files}")
    logger.info(f"[CN] Saved copy number in {cohort_name}.cn.tsv")
    logger.info(mergeCN(cn_files, cohort_name + ".cn.tsv"))

    # Allele Typing
    if not args.step_skip_typing:
        allele_files = alleleTyping(processed_bam, cn_files, method=args.allele_method)
        logger.debug(f"[Allele] Allele typing resuslt: {allele_files}")
        logger.info(f"[Allele] Saved in {cohort_name}.allele.tsv")
        mergeAllele(allele_files, cohort_name + ".allele.tsv")
    logger.info("[Main] Success")
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
