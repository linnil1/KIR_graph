from glob import glob
from pathlib import Path
from collections.abc import Iterable
from typing import Iterator, DefaultDict, Any
from collections import defaultdict
from dataclasses import dataclass
from itertools import chain
import os
import json

import pysam
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, dcc, html
from plotly.subplots import make_subplots

from namepipe import NamePath
from graphkir.plot import showPlot
from graphkir.utils import runShell
from graphkir.kir_cn import bam2Depth, readSamtoolsDepth


@dataclass(order=True)
class PositionInfo:
    ref: str
    pos: int
    depth: int


RegionPosition = tuple[str, list[PositionInfo]]


def extractCurrentAndPreviousRealPosition(
    bam_file: str,
) -> Iterator[tuple[str, int, str, str]]:
    """
    We extract the previous mapping info in read id
    with the format defined in bam2fastq function.

    The bamfile MUST derive from bam2fastq steps.

    yield: reference_mapped_on, pos_mapped_on, previous_reference, previous_position
    """
    for record in pysam.AlignmentFile(bam_file, "rb").fetch(until_eof=True):
        if not record.is_proper_pair:
            continue
        id = str(record.query_name)
        prev_mapped = id.split("|")[1:5]
        ref_name = str(record.reference_name)
        ref_pos = record.reference_start
        if record.is_read1:
            yield (ref_name, ref_pos, prev_mapped[0], prev_mapped[1])
        else:
            yield (ref_name, ref_pos, prev_mapped[2], prev_mapped[3])


def extractCurrentAndPreviousSimuPosition(
    bam_file: str,
) -> Iterator[tuple[str, int, str, int]]:
    """
    We extract the previous mapping info in read id.
    The read id will save the allele where the read generated from.
    (by art_illumina).

    previous_position is always 0 becuase the mapping position is stored in samfile not in id.

    yield: previous_reference, previous_position, reference_mapped_on, pos_mapped_on
    """
    for record in pysam.AlignmentFile(bam_file, "rb").fetch(until_eof=True):
        if not record.is_proper_pair:
            continue
        id = str(record.query_name)
        prev_ref = str(record.query_name).split("*")[0]
        ref_name = str(record.reference_name)
        ref_pos = record.reference_start
        yield (prev_ref, 0, ref_name, ref_pos)


def removePrevUnmapped(df: pd.DataFrame) -> pd.DataFrame:
    """Print and Remove previous unmapped reads"""
    print("Total")
    print(len(df))
    print("Unmapped")
    print(df[df["prev_pos"] == "*"])
    df = df[df["prev_pos"] != "*"]
    return df


def positionsSummary(
    positions: list[PositionInfo], depth_multiply: float = 1
) -> dict[str, Any]:
    """Summarize a list of depths and positions"""
    left = min(map(lambda i: i.pos, positions))
    right = max(map(lambda i: i.pos, positions))
    total_depth = sum(map(lambda i: i.depth, positions))
    length = right - left
    return {
        "ref": positions[0].ref,
        "left": left,
        "right": right,
        "length": length,
        "avg_depth": total_depth / length * depth_multiply,
    }


def extractSignificantRegion(
    positions: Iterable[PositionInfo],
    min_depth: float = 20,
    min_length: int = 20,
    max_gap: int = 5000,
) -> Iterator[list[PositionInfo]]:
    """ Find the continuous positions if the read depth is significant enough.  """
    ref_pos_dict = defaultdict(list)
    for pos_item in positions:
        ref_pos_dict[pos_item.ref].append(pos_item)

    for ref in ref_pos_dict:
        pos_depths = sorted(ref_pos_dict[ref])
        acc_positions: list[PositionInfo] = []
        total_depth = 0
        for pos_item in pos_depths:
            # extend
            if not acc_positions or pos_item.pos - acc_positions[-1].pos < max_gap:
                acc_positions.append(pos_item)
                total_depth += pos_item.depth
                continue
            # discard insignificant
            left = min(map(lambda i: i.pos, acc_positions))
            right = max(map(lambda i: i.pos, acc_positions))
            if (
                left == right
                or right - left < min_length
                or total_depth / (right - left) < min_depth
            ):
                acc_positions = []
                total_depth = 0
                continue
            # return
            yield acc_positions
            acc_positions = []
            total_depth = 0

        # the last one
        if len(acc_positions):
            left = min(map(lambda i: i.pos, acc_positions))
            right = max(map(lambda i: i.pos, acc_positions))
            if (
                left == right
                or right - left < min_length
                or total_depth / (right - left) < min_depth
            ):
                acc_positions = []
                total_depth = 0
                continue
            # return
            yield acc_positions


def extractPerReadMapping(cohort: str) -> Iterator[pd.DataFrame]:
    """Extract mapping infomation in bam to dataframe"""
    for input_name in NamePath(cohort).get_input_names():
        reads = extractCurrentAndPreviousRealPosition(input_name + ".bam")
        df = pd.DataFrame(reads, columns=["ref", "pos", "prev_ref", "prev_pos"])  # type: ignore
        df = removePrevUnmapped(df)
        df["prev_pos"] = df["prev_pos"].astype(int)
        yield df


def extractSignificantMappingRegion(
    df_reads: Iterable[pd.DataFrame],
) -> Iterator[RegionPosition]:
    """ Group mapping infomation by gene and yield the significant regions from them """
    # remove these three line to run individually
    df_reads_list = list(df_reads)
    df_reads = [pd.concat(df_reads_list)]  # merge all samples in cohort
    N = len(df_reads_list)
    for df in df_reads:
        # N = 1
        for gene in sorted(set(df["ref"])):
            df_gene = df[df["ref"] == gene]
            reads = map(
                lambda i: PositionInfo(ref=i.prev_ref, pos=i.prev_pos, depth=1),
                df_gene.itertuples(),
            )
            region_pos_item = extractSignificantRegion(reads, min_depth=N * 5 / 150)
            for pos_items in region_pos_item:
                yield gene.split("*")[0], pos_items


def printSignigicantRegion(
    regions: Iterable[RegionPosition], depth_multiply: float = 1
) -> pd.DataFrame:
    """Print the summary of the significant region"""
    mapping_info: list[dict[str, Any]] = []
    for gene_mapped_on, positions in regions:
        mapping_info.append(
            {
                **positionsSummary(positions, depth_multiply=depth_multiply),
                "gene": gene_mapped_on,
            }
        )
    df = pd.DataFrame(mapping_info)
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # type: ignore
        print(df)
        if len(df[df["gene"] == "KIR3DL2"]) > 10:  # some reads always mapped on KIR3DL2
            print(df[df["gene"] != "KIR3DL2"])
    return df


def plotSignigicantRegion(regions: Iterable[RegionPosition]) -> list[go.Figure]:
    """Plot the significant region"""
    gene_regions = defaultdict(list)
    for gene_mapped_on, positions in regions:
        gene_regions[gene_mapped_on].append(positions)

    figs: list[go.Figure] = []
    for gene_mapped_on, gene_positions in gene_regions.items():
        figs_gene = []
        for positions in gene_positions:
            fig = go.Histogram(x=list(map(lambda i: i.pos, positions)))
            fig.xbins.size = 150
            figs_gene.append((fig, positions[0].ref))

        if not len(figs_gene):
            continue
        fig_merge = make_subplots(
            rows=1, cols=len(figs_gene), shared_yaxes=True, horizontal_spacing=0.00
        )
        for i, (fig, ref) in enumerate(figs_gene):
            fig_merge.add_trace(fig, row=1, col=i + 1)
            fig_merge.update_xaxes(title_text=ref, row=1, col=i + 1)
        fig_merge.update_layout(
            title=gene_mapped_on,
            showlegend=False,
        )
        figs.append(fig_merge)

    # plot on same row
    return figs


def getNCBIKirGene(ref: str = "hg19") -> pd.DataFrame:
    """Get KIR genes position from NCBI record"""
    # read NCBI data
    df_ncbi = pd.read_csv("ncbi_kir_search_result.tsv", sep="\t")
    df_ncbi = df_ncbi[df_ncbi["Symbol"].str.contains("KIR")]
    print(df_ncbi[["Symbol", "GeneID"]])

    # read NCBI gene data
    if not os.path.exists("ncbi_kir_gene.json"):
        data = {}
        import requests

        for row in df_ncbi[["Symbol", "GeneID"]].itertuples():
            id = row.GeneID
            gene = row.Symbol
            data[gene] = {
                "ncbi": requests.get(
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                    f"?db=gene&id={id}&rettype=json&retmode=json"
                ).json(),
                "id": id,
                "gene": gene,
            }
        json.dump(data, open("ncbi_kir_gene.json", "w"))
    else:
        data = json.load(open("ncbi_kir_gene.json"))

    # extract hg19 hg38 position
    locations = []
    for gene, gene_info in data.items():
        id = gene_info["id"]
        regions = gene_info["ncbi"]["result"][str(id)]["locationhist"]
        release = sorted(set(i["annotationrelease"] for i in regions))
        if ref == "hg38":
            release = [i for i in release if i.startswith("110")]  # hg38
        else:
            release = [i for i in release if i.startswith("105")]  # hg19
        if not len(release):
            continue
        latest_release = release[-1]
        location = [i for i in regions if i["annotationrelease"] == latest_release]
        for loc in location:
            loc["gene"] = gene
        locations.extend(location)

    # rename and print
    df_loc = pd.DataFrame(locations)
    df_loc = df_loc.rename(
        columns={
            "chraccver": "reference",
            "chrstart": "pos_left",
            "chrstop": "pos_right",
        }
    )
    df_loc = df_loc[["gene", "reference", "pos_left", "pos_right"]]
    df_loc = df_loc.sort_values("gene")
    print("NCBI KIR gene data")
    print(df_loc)
    print(df_loc[df_loc["reference"].str.startswith("NC_")])
    return df_loc


def getUCSCKirGene(ref: str = "hg19") -> pd.DataFrame:
    """Get UCSC annotation of KIR gene"""
    if ref == "hg19":
        if not os.path.exists("hg19.refGene.gtf.gz"):
            runShell(
                "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz"
            )
            runShell(
                'zcat hg19.refGene.gtf.gz | grep "\\"KIR" | grep -v "exon" | grep -v "KIRREL" > hg19.refGene.kir.gtf'
            )
        gff = "hg19.refGene.kir.gtf"
    elif ref == "hg38":
        if not os.path.exists("hg38.refGene.gtf.gz"):
            runShell(
                "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz"
            )
            runShell(
                'zcat hg38.refGene.gtf.gz | grep "\\"KIR" | grep -v "exon" | grep -v "KIRREL" > hg38.refGene.kir.gtf'
            )
        gff = "hg38.refGene.kir.gtf"

    # to dataframe
    df = pd.read_csv(
        gff,
        sep="\t",
        names=[
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ],
    )
    df["gene"] = (
        df["attribute"].str.split(";", expand=True)[0].str.split('"', expand=True)[1]
    )
    df_gene_regions = (
        df[["gene", "seqname", "start", "end"]].sort_values("gene").drop_duplicates()
    )
    df_gene_regions.columns = ["gene", "reference", "pos_left", "pos_right"]  # type: ignore
    print("UCSC KIR gene")
    print(df_gene_regions)
    return df_gene_regions


def extractReadDepth(cohort: str) -> Iterator[pd.DataFrame]:
    """Read depth of bam file from samtools depth"""
    for input_name in NamePath(cohort).get_input_names():
        output_name = input_name + ".depth"
        if not Path(output_name + ".tsv").exists():
            bam2Depth(input_name + ".bam", output_name + ".tsv", get_all=False)
        df = readSamtoolsDepth(output_name + ".tsv")
        yield df


def extractSignificantDepthRegion(
    df_depths: Iterable[pd.DataFrame],
) -> Iterator[RegionPosition]:
    """Find the significant region from depth file"""
    df_depths = [pd.concat(df_depths)]
    for df in df_depths:
        df = df[df["depth"] > 5]
        print(df)
        depths = map(
            lambda i: PositionInfo(ref=str(i.gene), pos=i.pos, depth=i.depth),
            df.itertuples(),
        )
        for positions in extractSignificantRegion(depths):
            yield "KIR", positions


def extractSimulateMapping(cohort: str) -> Iterator[pd.DataFrame]:
    """Transform simulated bamfile to dataframe"""
    for input_name in NamePath(cohort).get_input_names():
        reads = extractCurrentAndPreviousSimuPosition(input_name + ".bam")
        df = pd.DataFrame(reads, columns=["ref", "pos", "prev_ref", "prev_pos"])  # type: ignore
        yield df


def evaluateSimulationReadMappingOnGenome(cohort: str) -> list[go.Figure]:
    figs = []
    # SELECT regions FROM ALL
    dfs_depth = extractReadDepth(cohort)
    regions = list(extractSignificantDepthRegion(dfs_depth))
    printSignigicantRegion(regions)
    figs.extend(plotSignigicantRegion(regions))

    # SELECT regions FROM ALL GROUP BY GENE
    dfs_read = extractSimulateMapping(cohort)
    regions = list(extractSignificantMappingRegion(dfs_read))
    printSignigicantRegion(regions, 150)
    figs.extend(plotSignigicantRegion(regions))
    return figs


def evaluateRealMappingOnKIR(cohort: str) -> list[go.Figure]:
    figs = []
    # SELECT regions FROM ALL GROUP BY GENE
    dfs_read = extractPerReadMapping(cohort)
    regions = list(extractSignificantMappingRegion(dfs_read))
    printSignigicantRegion(regions, 150)
    figs.extend(plotSignigicantRegion(regions))
    return figs


if __name__ == "__main__":
    figs: list[go.Figure] = []
    # simulate dataset
    cohort = "data6/linnil1_syn_s44.{}.30x_s444.index5_hs37d5.bwa"
    # figs.extend(evaluateSimulationReadMappingOnGenome(cohort))

    # real dataset
    ref = "hg19"
    # data in dataset
    # getNCBIKirGene(ref)
    # getUCSCKirGene(ref)

    # twbb 14 samples hg19 bam2fastq data
    cohort = "data_real/hg19.twbb.{}.part_merge.annot_read.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim"
    # figs.extend(evaluateRealMappingOnKIR(cohort))

    # twbb 14 samples hs37d5 bam2fastq data
    cohort = "data_real/twbb.{}.index_hs37d5.bwa.part_merge.annot_read.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim"
    # figs.extend(evaluateRealMappingOnKIR(cohort))

    # hprc 28 samples hs37d5 bam2fastq data
    cohort = "data_real/hprc.{}.index_hs37d5.bwa.part_merge.annot_read.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim"
    # figs.extend(evaluateRealMappingOnKIR(cohort))

    # hprc 28 samples hs37d5 with very limited extracted reads
    cohort = "data_real/hprc.{}.index_hs37d5.bwa.part_strict.annot_read.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim"
    # figs.extend(evaluateRealMappingOnKIR(cohort))

    # plot all figure
    showPlot(figs)
