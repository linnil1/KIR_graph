import os
import json
from glob import glob
from typing import Iterator, Iterable
from itertools import chain

import pysam
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, dcc, html
from graphkir.plot import showPlot
from plotly.subplots import make_subplots


def extractCurrentAndPreviousPosition(
    bam_file: str,
) -> Iterator[tuple[str, int, str, str]]:
    """yield: reference_mapped_on, pos_mapped_on, previous_reference, previous_position"""
    for record in pysam.AlignmentFile(bam_file, "rb").fetch(until_eof=True):
        # print(record.reference_name)
        # print(record.reference_name)
        if not record.is_proper_pair:
            continue
        # print(str(record))
        id = str(record.query_name)
        prev_mapped = id.split("|")[1:5]
        ref_name = str(record.reference_name)
        ref_pos = record.reference_start
        if record.is_read1:
            yield (ref_name, ref_pos, prev_mapped[0], prev_mapped[1])
        else:
            yield (ref_name, ref_pos, prev_mapped[2], prev_mapped[3])


def removePrevUnmapped(df: pd.DataFrame) -> pd.DataFrame:
    """Print and Remove previous unmapped reads"""
    print("Total")
    print(df)
    print("Unmapped")
    print(df[df["prev_pos"] == "*"])
    df = df[df["prev_pos"] != "*"]
    return df


def extractSignificantRegion(
    df_gene: pd.DataFrame, significant_count: int = 20
) -> Iterator[tuple[str, list[int]]]:
    """
    Remove outlier and split positions if they are too far away

    Yield: reference_name, positions
    """
    max_span_size = 5000
    for ref in sorted(set(df_gene["prev_ref"])):
        df_gene_ref = df_gene[df_gene["prev_ref"] == ref]
        x = sorted(list(df_gene_ref["prev_pos"]))
        current: list[int] = []
        for i in x:
            # extend
            if not current or i - current[-1] < max_span_size:
                current.append(i)
                continue
            # discard insignificant
            if len(current) < significant_count:
                continue
            # return
            yield ref, current
            current = []
    # the last one
    if len(current) >= significant_count:
        yield ref, current


def plotAllHistogram(rps: Iterable[tuple[str, list[int]]]) -> go.Figure | None:
    """Plot all mapping positions into one histogram"""
    # generate histogram
    figs_gene = []
    for reference, positions in rps:
        fig = go.Histogram(x=positions)
        fig.xbins.size = 150
        figs_gene.append((fig, reference))

    # check
    if not len(figs_gene):
        return None

    # plot on same row
    fig_all = make_subplots(
        rows=1, cols=len(figs_gene), shared_yaxes=True, horizontal_spacing=0.00
    )
    for i, (fig, ref) in enumerate(figs_gene):
        fig_all.add_trace(fig, row=1, col=i + 1)
        fig_all.update_xaxes(title_text=ref, row=1, col=i + 1)
    return fig_all


def printPrevMappingSummary(df: pd.DataFrame):
    """ print the previous mapping positions with current mapping positions """
    significant_count = 30
    mapping_info = []
    for gene in sorted(set(df["ref"])):
        df_gene = df[df["ref"] == gene]
        for reference, positions in extractSignificantRegion(
            df_gene, significant_count=significant_count
        ):
            mapping_info.append(
                {
                    "gene": gene,
                    "reference": reference,
                    "pos_left": min(positions),
                    "pos_right": max(positions),
                    "pos_num": len(positions),
                }
            )
    df_mapping_info = pd.DataFrame(mapping_info)
    print("Our mapping data")
    print(df_mapping_info)
    print(df_mapping_info[df_mapping_info["reference"] == "chr19"])


def plotPrevMappingSummary(df: pd.DataFrame):
    """ Plot the previous mapping positions with current mapping positions """
    # plot part
    significant_count = 30
    figs = []
    for gene in sorted(set(df["ref"])):
        print(gene)
        df_gene = df[df["ref"] == gene]
        fig = plotAllHistogram(
            extractSignificantRegion(df_gene, significant_count=significant_count)
        )
        if fig is None:
            continue
        fig.update_layout(title=gene)
        figs.append(fig)
    return figs


def plotPrevMapping() -> list[go.Figure]:
    """ Plot previous mapping with current gene"""
    bam_files = glob(
        # "data_twbb/hg19.twbb.*.part_merge.annot_read.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim.bam"
        "data_twbb/hg38.twbb.*.part_merge.annot_read.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim.bam"
    )
    print(bam_files)
    # bam_files = ["data_twbb/hg19.twbb.NGS2_20150412B.part_merge.annot_read.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim.bam"]
    # bam_file = "data_twbb/hg19.twbb.NGS2_20150505B.part_merge.annot_read.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim.bam"
    # data = extractCurrentAndPreviousPosition(bam_file)
    datas = []
    for bam_file in bam_files:
        datas.append(extractCurrentAndPreviousPosition(bam_file))
    data = chain.from_iterable(datas)

    # data part
    df = pd.DataFrame(data, columns=["ref", "pos", "prev_ref", "prev_pos"])  # type: ignore
    df = removePrevUnmapped(df)
    df["prev_pos"] = df["prev_pos"].astype(int)

    # data
    printPrevMappingSummary(df)
    # plot
    return plotPrevMappingSummary(df)


def getNCBIKirGene(ref: str = "hg38"):
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


if __name__ == "__main__":
    getNCBIKirGene()
    # work with twbb1 bam2fastq data
    showPlot(plotPrevMapping())
