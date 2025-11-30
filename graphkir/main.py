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
from .wgs import downloadHg19, bwa, bwaIndex, extractDiploidCoverage, extractFromHg19, bam2fastq
from .hisat2 import hisatMap, extractVariantFromBam, readExons
from .kir_cn import predictSamplesCN, loadCN, filterDepth
from .samtools_utils import bam2Depth
from .kir_typing import selectKirTypingModel
from .utils import setThreads, getThreads, mergeCN, mergeAllele, logger
from .external_tools import setEngine
from .plot import plotCN, plotReadMappingStat, showPlot, plotGeneDepths, savePlot


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
        names: list[str], reads: list[tuple[str, str]], index_wgs: str, diploid_gene: str
) -> tuple[list[str], list[tuple[str, str]], list[str]]:
    new_names = []
    new_reads = []
    diploid_depths = []
    for name, (fq1, fq2) in zip(names, reads):
        # read mapping
        logger.info(f"[WGS] Run BWA mapping with index {index_wgs} on {name}")
        suffix = "." + index_wgs.replace(".", "_").replace("/", "_")
        # bwa(index_wgs, fq1, fq2, name + suffix, threads=getThreads())
        name += suffix

        # extract diploid coverage information
        if diploid_gene != '':
            diploid_depths.append(extractDiploidCoverage(name, diploid_gene))
        else:
            diploid_depths.append('')

        # extract KIR
        suffix = ".extract"
        # logger.info(f"[WGS] Extract chr19... from {name}")  # written in func
        extractFromHg19(name + ".bam", name + suffix, "hs37d5", threads=getThreads())
        name += suffix
        logger.info(f"[WGS] Extract read from {name}.bam")
        bam2fastq(name + ".bam", name, threads=getThreads())
        new_names.append(name)
        new_reads.append((name + ".read.1.fq.gz", name + ".read.2.fq.gz"))
    return new_names, new_reads, diploid_depths


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
        extractVariantFromBam(index_ref, name + ".bam", name + suffix, error_correction=False)
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
        called_alleles, warning_genes = t.typing(cn)
        logger.info(f"[Allele] {called_alleles} ({name})")
        name += suffix
        # t.save(name + ".json")
        df = pd.DataFrame(
            {
                "name": [name],
                "alleles": ["_".join(called_alleles)],
                "warnings": ["_".join(warning_genes)],
            }
        )
        df.to_csv(name + ".tsv", sep="\t", index=False)
        allele_files.append(name + ".tsv")

        logger.info(
            f"[Allele] All possible allele set in Allele typing saved in [name].possible.tsv"
        )
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
        description="Run Graph-KIR",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # General Options
    parser.add_argument(
        "--thread",
        default=1,
        help="Number of threads",
    )
    parser.add_argument(
        "--engine",
        default="local",
        choices=["podman", "docker", "singularity", "local"],
        help="Run the external package with podman, docker, or singularity. "
        "Alternatively, run it in a local environment, which requires manual installation.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=logging._nameToLevel.keys(),
        help="Set the log level",
    )

    # Input/Output Data
    parser.add_argument(
        "--r1",
        action="append",
        help="Paths to paired-end Read 1 FASTQ files (provide in order for multiple samples)",
    )
    parser.add_argument(
        "--r2",
        action="append",
        help="Paths to paired-end Read 2 FASTQ files (in order)",
    )
    parser.add_argument(
        "--input-csv",
        help="Path to a CSV file containing samples (columns: id, fastq, cnfile)",
    )
    parser.add_argument(
        "--output-folder",
        help="Output folder for saving data (default: same folder as input reads)",
    )
    parser.add_argument(
        "--output-cohort-name",
        help="Output prefix to save the result of the entire cohort (Default: {output-folder}/cohort)",
    )

    parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate plots for intermediate results (saved as {output-cohort-name}.xx.png)",
    )

    # Index and MSA
    parser.add_argument(
        "--ipd-version",
        default="2100",
        help="IPD*KIR version. Branch name in https://github.com/ANHIG/IPDKIR are all available.",
    )
    parser.add_argument(
        "--msa-type",
        default="ab_2dl1s1",
        choices=["merge", "split", "ab", "ab_2dl1s1"],
        help="Type of Multiple Sequence Alignment (MSA) setup for Graph-KIR. "
        "Options: merge (merge KIR genes into 1 gene), split (17 genes), "
        "ab (16 genes, merging 2DL5A and 2DL5B), ab_2dl1s1 (15 genes, merging 2DL1 and 2DS1).",
    )
    parser.add_argument(
        "--msa-no-exon-only-allele",
        action="store_true",
        help="Exclude exon-only alleles from the allele set. (requires index rebuilding)",
    )
    parser.add_argument(
        "--index-folder",
        default="index",
        help="The path to the index folder, which must include the HISAT2-indexed KIR. "
        "Optionally, the hs37d5 index can also be located in the same folder. "
        "Alternatively, you can specify the path to the hs37d5 index using `--index-wgs`.",
    )
    parser.add_argument(
        "--index-wgs",
        help="Path to a BWA-indexed hs37d5 index file (Default=index_folder/hs37d5.fa.gz)",
    )

    # Copy Number
    parser.add_argument(
        "--cn-diploid-gene",
        choices=['', "VDR", "RYR1", "EGFR"],
        default='',
        help="Specify a gene to provide diploid coverage information during the CN estimation step. "
        "(Default is '', which means the CN estimation module will operate without diploid information.) "
        "Options for reference diploid genes are: VDR, RYR1, and EGFR."
    )
    parser.add_argument(
        "--cn-exon",
        action="store_true",
        help="Select exon-only regions of genes for copy number (CN) prediction instead of all positions",
    )
    parser.add_argument(
        "--cn-cohort",
        action="store_true",
        help="Predict CN cohort-wise instead of sample-wise",
    )
    parser.add_argument(
        "--cn-select",
        default="p75",
        choices=["p75", "mean", "median"],
        help="Method for selecting depths of each gene; 75 percentile(p75), mean(mean), median(median).",
    )
    parser.add_argument(
        "--cn-algorithm",
        default="LCND",
        choices=["LCND", "KDE"],
        help="CN prediction model. Choose between LCND (Linear Copy Number Distribution, the method described in the paper) or KDE (Kernel Density Estimation).",
    )
    parser.add_argument(
        "--cn-dist-dev",
        default=0.08,
        help="Deviation of distrubionts in LCND (sometimes 0.06 works better)",
    )
    parser.add_argument(
        "--cn-3dl3-not-diploid",
        action="store_true",
        help="Do not assume KIR3DL3 as a diploid gene (for LCND)",
    )
    parser.add_argument(
        "--cn-provided",
        nargs="*",
        help="Provided CN predictions in TSV format (each file represents a sample in order)",
    )

    parser.add_argument(
        "--allele-strategy",
        default="full",
        choices=["full", "exonfirst", "report"],
        help="Choose the allele typing strategy: 'full' for maximum likelihood using all variants, "
        "'exonfirst' for typing exon variants before full variants, "
        "or 'report' for typing via abundance using an EM-algorithm similar to HISAT2-genotype.",
    )

    # Skip Processing Steps
    parser.add_argument(
        "--step-skip-extraction",
        action="store_true",
        help="Skip filtering KIR reads from the (WGS) fastq files",
    )
    parser.add_argument(
        "--step-skip-typing",
        action="store_true",
        help="Skip allele-typing",
    )
    return parser


def main(args: argparse.Namespace) -> None:
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
        if not args.r1:
            raise ValueError("At least one paired-end read 1 FASTQ file must be provided")
        if len(args.r1) != len(args.r2):
            raise ValueError("The number of paired-end read 1 and read 2 FASTQ files must be equal")
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
        raise ValueError("No samples found in input")
    if len(cn_files) != len(names):
        raise ValueError("Mismatch between number of samples and copy number files")
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
        diploid_gene = args.cn_diploid_gene
        if args.cn_cohort:
            diploid_gene = ""
        names, reads, diploid_depths = runWGS(names, reads, index_wgs, diploid_gene)
    else:
        diploid_depths = ['' for _ in range(len(names))]



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
        "base_dev": float(args.cn_dist_dev),
        "start_base": 2,
    }
    # Copy Number determination
    if all(cn_files):
        pass
    elif not args.cn_cohort:
        # per sample
        for i, depth_file in enumerate(depth_files):
            if cn_files[i]:  # CN file exists, skip
                continue
            suffix = f".{args.cn_select}.{args.cn_algorithm}"
            name = str(Path(depth_file).with_suffix(suffix))
            diploid_depth = diploid_depths[i]
            logger.info(f"[CN] Copy number estimation per sample ({name})")
            predictSamplesCN(
                [depth_file],
                [name + ".tsv"],
                diploid_depth,
                cluster_method=args.cn_algorithm,
                cluster_method_kwargs=cluster_method_kwargs,
                assume_3DL3_diploid=not args.cn_3dl3_not_diploid,
                save_cn_model_path=name + ".json",
                select_mode=args.cn_select,
            )
            cn_files[i] = name + ".tsv"
            figs.extend(plotGeneDepths(depth_file))
            figs.extend(plotCN(name + ".json"))
    else:
        # by cohort
        suffix = f".{args.cn_select}.cohort.{args.cn_algorithm}"
        cn_cohort_name = cohort_name + suffix
        cn_files = [
            str(Path(path).with_suffix(suffix + ".tsv")) for path in depth_files
        ]
        logger.info(f"[CN] Copy number estimation by cohort ({cn_cohort_name})")
        predictSamplesCN(
            [depth_files[i] for i, cnf in enumerate(cn_files) if cnf],
            [cnf for i, cnf in enumerate(cn_files) if cnf],
            cluster_method=args.cn_algorithm,
            cluster_method_kwargs=cluster_method_kwargs,
            save_cn_model_path=cn_cohort_name + ".json",
            select_mode=args.cn_select,
        )
        # TODO: case of per gene CN prediction is plot from this json file
        figs.extend(plotCN(cn_cohort_name + ".json"))
    logger.debug(f"[CN] Copy number files: {cn_files}")
    logger.info(f"[CN] Saved copy number in {cohort_name}.cn.tsv")
    logger.info(mergeCN(cn_files, cohort_name + ".cn.tsv"))

    # Allele Typing
    if not args.step_skip_typing:
        allele_files = alleleTyping(
            processed_bam, cn_files, method=args.allele_strategy
        )
        logger.debug(f"[Allele] Allele typing resuslt: {allele_files}")
        logger.info(f"[Allele] Saved in {cohort_name}.allele.tsv")
        mergeAllele(allele_files, cohort_name + ".allele.tsv")
    logger.info("[Main] Success")

    if args.plot:
        savePlot(cohort_name + ".plot.html", figs)
        showPlot(figs)


def entrypoint() -> None:
    """Command line entry function"""
    parser = createParser()
    args = parser.parse_args()
    main(args)


if __name__ == "__main__":
    entrypoint()
