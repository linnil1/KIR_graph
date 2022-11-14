import os
import json
import subprocess
from pprint import pprint
from typing import Iterator, Any, TypedDict, Callable, Iterable
from itertools import chain
from dataclasses import dataclass
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor

from dash import Dash, dcc, html, Input, Output
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

from namepipe import NamePath
from graphkir.kir_cn import readSamtoolsDepth
from graphkir.hisat2 import loadReadsAndVariantsData
from graphkir.utils import samtobam, getGeneName, getAlleleField, runShell
from kg_utils import runDocker


@dataclass
class ReadRecord:
    id: str
    flag: int
    ref: str


class BamFile(TypedDict, total=False):
    id: int
    method: str
    file: str
    gene_compare_type: str
    records: dict[str, list[ReadRecord]]
    total: int


DfDict = dict[str, Any]


def extractAlleleName(seq: str) -> str:
    """ KIR2DL1*0010101-1-2123 83 ...  ->   KIR2DL1*0010101 """
    return seq.split("-")[0]


def samtools(cmd: str, name: str) -> Iterator[str]:
    """ Run samtools subcommand and read result """
    proc = runDocker("samtools", f"samtools {cmd} -@4 {name}", capture_output=True)
    return filter(None, proc.stdout.split("\n"))


def getStat(filename: str) -> dict[str, int]:
    """ Read total, secondary and properly paired from samtools stat """
    num_pair = 0
    num_total = 0
    num_second = 0
    for i in samtools("flagstat", filename):
        # 122050 + 0 properly paired (100.00% : N/A)
        # 122050 + 0 in total (QC-passed reads + QC-failed reads)
        # 0 + 0 secondary
        if "in total" in i:
            num_total = int(i.split()[0])
        if "secondary" in i:
            num_second = int(i.split()[0])
        if "properly paired" in i:
            num_pair = int(i.split()[0])
    return {
        "total": num_total,
        "secd": num_second,
        "pair": num_pair,
    }


def plotStat(df: pd.DataFrame) -> list[go.Figure]:
    """ Plot getStat info """
    ans_df = df[(df["method"] == "answer")]
    ans = dict(zip(ans_df["id"], ans_df["total"]))
    df["ans_total"] = df.apply(lambda i: ans[i.id], axis=1)
    df["pair_perc"] = df["pair"] / df["ans_total"]
    df["secd_perc"] = df["secd"] / df["ans_total"]
    display_columns = sorted(set(df["method"]))
    # display_columns = ["linear", "linear_ab", "linear_ab_2dl1s1", "linear_full",
    #                    "hisat", "hisat_ab", "hisat_ab_2dl1s1", "hisat_merge"]
    df = df[df["method"].isin(display_columns)]
    fig0 = px.box(df, x="method",  y="pair_perc", title="Proper Paired Ratio",
                  labels={"pair_perc": "Primary paired reads / Total reads"},
                  category_orders={'method': display_columns})
    fig1 = px.box(df, x="method",  y="secd_perc", title="Secondary Ratio",
                  labels={"secd_perc": "Secondary reads / Total reads"},
                  category_orders={'method': display_columns})
    return [fig0, fig1]


def plotBamStat() -> list[go.Figure]:
    # sample_index = "data/linnil1_syn_wide.test10"
    answer = "linnil1_syn/linnil1_syn_s2022.{}.30x_s1031.read..sam"
    prefix = "data6/linnil1_syn_s2022.{}.30x_s1031"
    filename = "tmp"
    # filename = f"{sample_index}.bam_stat_vg"

    cohort = [
        {"method": "answer",           "name": f"{answer}"},
        {"method": "linear",           "name": f"{prefix}.index5_kir_2100_withexon_split.leftalign.backbone.bowtie2.bam"},
        {"method": "bwa",              "name": f"{prefix}.index5_kir_2100_withexon_split.leftalign.backbone.bwa.bam"},
        {"method": "hisat",            "name": f"{prefix}.index5_kir_2100_withexon_split.leftalign.mut01.graph.bam"},
        {"method": "hisat_ab",         "name": f"{prefix}.index5_kir_2100_withexon_ab.leftalign.mut01.graph.bam"},
        {"method": "hisat_ab_2dl1s1",  "name": f"{prefix}.index5_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.bam"},
        # {"method": "linear_ab",                    "file": f"{name}.index_kir_2100_raw_cons.bowtie2.bam"},
        # {"method": "linear_ab_2dl1s1",             "file": f"{name}.index_kir_2100_2dl1s1_cons.bowtie2.bam"},
        # {"method": "linear_full",                  "file": f"{name}.index_kir_2100_raw_full.bowtie2.bam"},
        # {"method": "bwa_ab",                       "file": f"{name}.index_kir_2100_raw_cons_bwa.bwa.bam"},
        # {"method": "bwa_ab_2dl1s1",                "file": f"{name}.index_kir_2100_2dl1s1_cons_bwa.bwa.bam"},
        # {"method": "bwa_merge",                    "file": f"{name}.index_kir_2100_merge_cons_bwa.bwa.bam"},
        # {"method": "hisat_merge",                  "file": f"{name}.index_kir_2100_merge.mut01.bam"},
        # {"method": "hisat_merge_type",             "file": f"{name}.index_kir_2100_merge.mut01.hisatgenotype.errcorr.bam"},
        # {"method": "hisat_merge_type_nomulti",     "file": f"{name}.index_kir_2100_merge.mut01.hisatgenotype.errcorr.no_multi.bam"},
        # {"method": "hisat_ab_type",                "file": f"{name}.index_kir_2100_raw.mut01.hisatgenotype.errcorr.bam"},
        # {"method": "hisat_ab_type_nomulti",        "file": f"{name}.index_kir_2100_raw.mut01.hisatgenotype.errcorr.no_multi.bam"},
        # {"method": "hisat_ab_2dl1s1_type",         "file": f"{name}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.bam"},
        # {"method": "hisat_ab_2dl1s1_type_nomulti", "file": f"{name}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.no_multi.bam"},
        # {"method": "hisat_type",                   "file": f"{name}.index_kir_2100_ab.mut01.hisatgenotype.errcorr.bam"},
        # {"method": "hisat_type_nomulti",           "file": f"{name}.index_kir_2100_ab.mut01.hisatgenotype.errcorr.no_multi.bam"},
        # {"method": "hisat_271-ab-2dl1s1",          "file": f"{name}.index_kir_271_2dl1s1.mut01.bam"},
        # {"method": "hisat_290-ab-2dl1s1",          "file": f"{name}.index_kir_290_2dl1s1.mut01.bam"},
        # {"method": "vg_ab",                        "file": f"data3/{answer}.{id:02d}.data3_kir_2100_raw.vgindex.bam"},
        # {"method": "vg_merge",                     "file": f"data3/{answer}.{id:02d}.data3_kir_2100_merge.vgindex.bam"},
        # {"method": "ping_allele_setup",            "file": f"data3/ping_{answer}.result/allele_setup_files/{answer}.{id:02d}.read..bam"},
    ]

    if not os.path.exists(filename + ".json"):
        # Figure: Proper Paired Perecetage
        id = 0
        data: list[DfDict] = []
        for info in cohort:
            for name in NamePath(info["name"]).get_input_names()[:10]:
                data.append({
                    "method": info["method"],
                    "file": name,
                    "id": name.template_args[0],
                })
        pprint(data)

        # for i in data:
        #     i.update(getStat(i["file"]))
        with ProcessPoolExecutor() as executor:
            for d, result in zip(data, executor.map(getStat, map(lambda i: i["file"], data))):
                d.update(result)
        # json.dump(data, open(filename + ".json", "w"))
    else:
        data = json.load(open(filename + ".json"))

    df = pd.DataFrame(data)
    runShell("stty echo opost")
    runShell("stty echo opost")
    print(df)
    return plotStat(df)


def getEachReadMappedOn(filename: str) -> dict[str, list[ReadRecord]]:
    """
    Read bam file, and return where the read is mapped on
    in the form of dict[read-id, list[(flag, reference)]]
    """
    data = defaultdict(list)
    for i in samtools("view", filename):
        fields = i.split("\t")
        if len(fields) < 3:
            print(i)
            continue
        data[fields[0]].append(ReadRecord(
            id=fields[0],
            # split by "-" is for reading answer"s sam file
            flag=int(fields[1]),
            ref=fields[2].split("-")[0],
        ))
    return data


def customMissingCalc(total: dict[str, int],
                      reads: dict[str, list[ReadRecord]],
                      **kwargs: Any) -> Iterator[DfDict]:
    """
    Calculate the number of read cannot mapped

    Args:
      total: Nubmer(value) of reads belong to the reference(key)
      reads: The gene mapped on (value) of for each read(key)
    """
    data = {gene: {"total": num, "count": 0, "miss": 0} for gene, num in total.items()}

    print(sorted(reads.keys()))
    for read_name, mapping_info in reads.items():
        data[getGeneName(read_name)]["count"] += 1
        # don't include secondary, but it's find
        # because secondary occur only when multiple mapping.
        if len(list(filter(lambda i: i.ref != "*" and not (i.flag & 256), mapping_info))) < 2:
            data[getGeneName(read_name)]["miss"] += 1

    for gene, d in data.items():
        d["miss"] += d["total"] - d["count"]  # maybe removed by mapper
        yield {
            "gene": gene,
            "total": d["total"],
            "count": d["count"],
            "miss": d["miss"],
            "miss_perc": d["miss"] / d["total"],
        }


def customMissingPlot(df: pd.DataFrame) -> list[go.Figure]:
    gene_order = {"gene": sorted(set(df["gene"]))}
    return [
        px.box(df, x="gene", y="miss_perc", color="method", category_orders=gene_order)
          .update_layout(title="Missed reads belong to (Proportion)",
                         yaxis_title="Missed_reads / total_read_counts in that gene",
                         yaxis_tickformat="e"),
        px.box(df, x="method", y="miss_perc", color="method", category_orders=gene_order)
          .update_layout(title="Missed reads Per Method (Proportion)",
                         yaxis_title="Missed_reads / total_read_counts in that gene",
                         yaxis_tickformat="e"),
    ]


def customProperMappedCalc(total: dict[str, int],
                           reads: dict[str, list[ReadRecord]],
                           **kwargs: Any) -> Iterator[DfDict]:
    """ Any of pair-reads are correctly mapped on reference """
    data = {gene: {"total": num, "count": 0, "pair": 0} for gene, num in total.items()}

    for read_name, mapping_info in reads.items():
        data[getGeneName(read_name)]["count"] += 1
        if len(list(filter(lambda i: (i.flag & 2), mapping_info))) >= 2:
            data[getGeneName(read_name)]["pair"] += 1

    for gene, d in data.items():
        yield {
            "gene": gene,
            "total": d["total"],
            "count": d["count"],
            "pair_num": d["pair"],
            "pair_perc": d["pair"] / d["total"],
        }


def customProperMappedPlot(df: pd.DataFrame) -> list[go.Figure]:
    gene_order = {"gene": sorted(set(df["gene"]))}
    return [
        px.box(df, x="gene", y="pair_perc", color="method", category_orders=gene_order)
          .update_layout(title="Proper mapped",
                         yaxis_title="Correct read / total read (%)"),
        px.box(df, x="method", y="pair_perc", category_orders=gene_order)
          .update_layout(title="Proper mapped",
                         yaxis_title="Pair-mapped read / total read (%)")
    ]


def customSecondaryCalc(total: dict[str, int],
                        reads: dict[str, list[ReadRecord]],
                        **kwargs: Any) -> Iterator[DfDict]:
    """ The number of secondary reads """
    data = {gene: {"total": num, "count": 0, "secd": 0} for gene, num in total.items()}
    for read_name, mapping_info in reads.items():
        data[getGeneName(read_name)]["count"] += 1
        data[getGeneName(read_name)]["secd"] += len(list(filter(lambda i: (i.flag & 256), mapping_info)))

    for gene, d in data.items():
        yield {
            "gene": gene,
            "total": d["total"],
            "count": d["count"],
            "secd_num": d["secd"],
            "secd_perc": d["secd"] / d["count"],
        }


def customSecondaryPlot(df: pd.DataFrame) -> list[go.Figure]:
    gene_order = {"gene": sorted(set(df["gene"]))}
    return [
        px.box(df, x="gene", y="secd_perc", color="method", category_orders=gene_order)
          .update_layout(title="Secondary count",
                         yaxis_title="secondary read / total read (%)"),
        px.box(df, x="method", y="secd_perc", category_orders=gene_order)
          .update_layout(title="Secondary Count",
                         yaxis_title="secondary read / total read (%)")
    ]


def customGeneAccCalc(total: dict[str, int],
                      reads: dict[str, list[ReadRecord]],
                      getRenamedGeneName: Callable[[str], str],
                      **kwargs: Any) -> Iterator[DfDict]:
    """ The pair reads are mapped on reference """
    data = {gene: {"total": num, "count": 0, "unique": 0, "primary": 0, "secondary": 0} for gene, num in total.items()}

    for read_name, mapping_info in reads.items():
        data[getGeneName(read_name)]["count"] += 1
        # remove non-pair
        mapping_info = list(filter(lambda i: (i.flag & 2), mapping_info))

        primary = list(filter(lambda i: i.ref != "*" and not (i.flag & (256 | 2048)), mapping_info))
        # assert primary is 2
        if len(primary) not in [0, 2]:
            print(read_name, mapping_info)
            assert False
        if primary and getRenamedGeneName(primary[0].ref) == getRenamedGeneName(read_name):
            data[getGeneName(read_name)]["primary"] += 1

        second = list(filter(lambda i: i.ref != "*", mapping_info))
        if second and any(getRenamedGeneName(read.ref) == getRenamedGeneName(read_name) for read in second):
            data[getGeneName(read_name)]["secondary"] += 1

        # unique mapping
        unique_pair = list(filter(lambda i: i.ref != "*", mapping_info))
        if len(unique_pair) == 2 and getRenamedGeneName(primary[0].ref) == getRenamedGeneName(read_name):
            data[getGeneName(read_name)]["unique"] += 1

    for gene, d in data.items():
        yield {
            "gene": gene,
            "total": d["total"],
            "correct": d["primary"],
            "acc": d["primary"] / d["count"],
            "type": "primary",
        }
        yield {
            "gene": gene,
            "total": d["total"],
            "correct": d["secondary"],
            "acc": d["secondary"] / d["count"],
            "type": "primary+secondary",
        }
        yield {
            "gene": gene,
            "total": d["total"],
            "correct": d["unique"],
            "acc": d["unique"] / d["count"],
            "type": "unique",
        }


def customGeneAccPlot(df: pd.DataFrame) -> list[go.Figure]:
    gene_order = {"gene": sorted(set(df["gene"]))}
    group, color = reColor(df["method"], df["type"])
    return [
        px.box(df, x="gene", y="acc",
               color=group,
               color_discrete_map=color,
               category_orders=gene_order,
           )
           .update_layout(
               title="Gene-level accuracy(Primary)",
               yaxis_title="Correct read / total read (%)",
           ),
        px.box(df, x="method", y="acc", color="type", category_orders=gene_order)
          .update_layout(title="Gene-level accuracy Per Method (Proportion)",
                         yaxis_title="Correct read / total read (%)"),
    ]


def reColor(arr1: Iterable[str], arr2: Iterable[str]) -> tuple[list[str], dict[str, str]]:
    """
    Assign color for arr1_i x arr2_j

    Returns:
      group_name:  list of arr1_i x arr2_i
      name_to_color: key=arr1_i-arr2_j value=color
    """
    group_name = list(map(lambda i: i[0] + " - " + i[1], zip(arr1, arr2)))
    colors = (px.colors.qualitative.Dark2, px.colors.qualitative.Set2, px.colors.qualitative.Pastel2)

    group_color = {}
    for i, m in enumerate(sorted(set(arr1))):
        for j, t in enumerate(sorted(set(arr2))):
            group_color[m + " - " + t] = colors[j % len(colors[0])][i % len(colors[1])]
    return group_name, group_color


def findAnswerSample(data: Iterable[BamFile], id: int) -> Any:
    """ Find the sample dict with same id """
    for sample in data:
        if sample["method"] == "answer" and sample["id"] == id:
            return sample
    raise ValueError


def renameAb2dl1s1(name: str) -> str:
    name = getGeneName(name)
    return {
        "KIR2DL5A": "KIR2DL5",
        "KIR2DL5B": "KIR2DL5",
        "KIR2DL1": "KIR2DL1S1",
        "KIR2DS1": "KIR2DL1S1",
    }.get(name, name)


def applyStatFuncToBam(stat_func: Callable[..., Iterator[DfDict]],
                       data: list[BamFile]) -> pd.DataFrame:
    """ Calculate the stat by stat_func for each sample(data) """
    custom_stats: list[DfDict] = []
    for sample in data:
        if sample["method"] == "answer":
            continue
        print(sample["method"], sample["id"])
        records_ans = findAnswerSample(data, sample["id"])["records"]
        records = sample["records"]
        gene_ans_total = Counter(map(lambda name: getGeneName(name), records_ans.keys()))

        if sample["gene_compare_type"] == "":
            getRenamedGeneName: Callable[[str], str] = getGeneName
        elif sample["gene_compare_type"] == "ab":  # KIR2DL5A and KIR2DL5B
            getRenamedGeneName = lambda i: i[:7]
        elif sample["gene_compare_type"] == "ab_2dl1s1":
            getRenamedGeneName = renameAb2dl1s1

        stats = stat_func(gene_ans_total, records, getRenamedGeneName=getRenamedGeneName)
        custom_stats.extend([{
            "method": sample["method"],
            "id": sample["id"],
            **stat,
        } for stat in stats])

    df = pd.DataFrame(custom_stats)
    print(df)
    return df


def plotGenewiseMapping() -> list[go.Figure]:
    # sample_index = "data/linnil1_syn_wide.test10"
    answer = "linnil1_syn/linnil1_syn_s2022.{}.30x_s1031.read..sam"
    prefix = "data6/linnil1_syn_s2022.{}.30x_s1031"
    answer = "linnil1_syn/linnil1_syn_s44.{}.30x_s444.read..sam"
    prefix = "data6/linnil1_syn_s44.{}.30x_s444"
    filename = "tmp"
    # filename = f"{sample_index}.bam_stat_vg"

    cohort = [
        {"method": "answer",          "compare_gene": "",          "name": f"{answer}"},
        # {"method": "linear",          "compare_gene": "",          "name": f"{prefix}.index5_kir_2100_withexon_split.leftalign.backbone.bowtie2.bam"},
        # {"method": "bwa",             "compare_gene": "",          "name": f"{prefix}.index5_kir_2100_withexon_split.leftalign.backbone.bwa.bam"},
        {"method": "hisat",           "compare_gene": "",          "name": f"{prefix}.index5_kir_2100_withexon_split.leftalign.mut01.graph.bam"},
        {"method": "hisat_ab",        "compare_gene": "ab",        "name": f"{prefix}.index5_kir_2100_withexon_ab.leftalign.mut01.graph.bam"},
        {"method": "hisat_ab_2dl1s1", "compare_gene": "ab_2dl1s1", "name": f"{prefix}.index5_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.bam"},
        # {"method": "linear_ab",                    "file": f"{name}.index_kir_2100_raw_cons.bowtie2.bam"},
        # {"method": "linear_ab_2dl1s1",             "file": f"{name}.index_kir_2100_2dl1s1_cons.bowtie2.bam"},
        # {"method": "linear_full",                  "file": f"{name}.index_kir_2100_raw_full.bowtie2.bam"},
        # {"method": "bwa_ab",                       "file": f"{name}.index_kir_2100_raw_cons_bwa.bwa.bam"},
        # {"method": "bwa_ab_2dl1s1",                "file": f"{name}.index_kir_2100_2dl1s1_cons_bwa.bwa.bam"},
        # {"method": "bwa_merge",                    "file": f"{name}.index_kir_2100_merge_cons_bwa.bwa.bam"},
        # {"method": "hisat_merge",                  "file": f"{name}.index_kir_2100_merge.mut01.bam"},
        # {"method": "hisat_merge_type",             "file": f"{name}.index_kir_2100_merge.mut01.hisatgenotype.errcorr.bam"},
        # {"method": "hisat_merge_type_nomulti",     "file": f"{name}.index_kir_2100_merge.mut01.hisatgenotype.errcorr.no_multi.bam"},
        # {"method": "hisat_ab_type",                "file": f"{name}.index_kir_2100_raw.mut01.hisatgenotype.errcorr.bam"},
        # {"method": "hisat_ab_type_nomulti",        "file": f"{name}.index_kir_2100_raw.mut01.hisatgenotype.errcorr.no_multi.bam"},
        # {"method": "hisat_ab_2dl1s1_type",         "file": f"{name}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.bam"},
        # {"method": "hisat_ab_2dl1s1_type_nomulti", "file": f"{name}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.no_multi.bam"},
        # {"method": "hisat_type",                   "file": f"{name}.index_kir_2100_ab.mut01.hisatgenotype.errcorr.bam"},
        # {"method": "hisat_type_nomulti",           "file": f"{name}.index_kir_2100_ab.mut01.hisatgenotype.errcorr.no_multi.bam"},
        # {"method": "hisat_271-ab-2dl1s1",          "file": f"{name}.index_kir_271_2dl1s1.mut01.bam"},
        # {"method": "hisat_290-ab-2dl1s1",          "file": f"{name}.index_kir_290_2dl1s1.mut01.bam"},
        # {"method": "vg_ab",                        "file": f"data3/{answer}.{id:02d}.data3_kir_2100_raw.vgindex.bam"},
        # {"method": "vg_merge",                     "file": f"data3/{answer}.{id:02d}.data3_kir_2100_merge.vgindex.bam"},
        # {"method": "ping_allele_setup",            "file": f"data3/ping_{answer}.result/allele_setup_files/{answer}.{id:02d}.read..bam"},
    ]

    if not os.path.exists(filename + ".json"):
        # Figure: Proper Paired Perecetage
        id = 0
        data: list[DfDict] = []
        for info in cohort:
            for name in NamePath(info["name"]).get_input_names()[:1]:
                data.append({
                    "gene_compare_type": info["compare_gene"],
                    "method": info["method"],
                    "file": name,
                    "id": name.template_args[0],
                })
        pprint(data)

        with ProcessPoolExecutor() as executor:
            bam_files = map(lambda i: i["file"], data)
            for d, result in zip(data, executor.map(getEachReadMappedOn, bam_files)):
                d["records"] = result  # reocrds by parseing samtools view
                d["total"] = len(result)

        # print(data)
        # json.dump(data, open(filename + ".json", "w"))
    else:
        print("Reading")
        data = json.load(open(filename + ".json"))

    figs = []
    # missing
    df = applyStatFuncToBam(customMissingCalc, data)
    figs.extend(customMissingPlot(df))
    # Any proper reads
    df = applyStatFuncToBam(customProperMappedCalc, data)
    figs.extend(customProperMappedPlot(df))
    # secondary
    df = applyStatFuncToBam(customSecondaryCalc, data)
    figs.extend(customSecondaryPlot(df))
    # gene-acc
    df = applyStatFuncToBam(customGeneAccCalc, data)
    figs.extend(customGeneAccPlot(df))
    return figs


def isGeneCorrect(a: str, b: str) -> bool:
    """
    Compare two gene are the same

    This function handle the special case of KIR2DL5 and KIR2DL1S1
    """
    if a == "KIR2DL5":
        # b = "KIR2DL5A" or "KIR2DL5B"
        return a in b
    elif a == "KIR2DL1S1":
        return b in ["KIR2DL1", "KIR2DS1"]
    return a == b


def calcFromToPerSample(filename: str) -> pd.DataFrame:
    """ extract the gene came from and the gene the reads mapped on per reads """
    reads = loadReadsAndVariantsData(filename)["reads"]
    data_from_to = []
    for read in reads:
        gene = getGeneName(read.backbone)
        ans_gene = getGeneName(read.l_sam)
        data_from_to.append({
            "from_allele": extractAlleleName(read.l_sam),
            "from": ans_gene,
            "to": gene,
            "correct": isGeneCorrect(gene, ans_gene),
            "is_multi": read.multiple > 1,
            "multi": read.multiple,
        })
    return pd.DataFrame(data_from_to)


def plotFromToBar(df: pd.DataFrame, x: str, y: str) -> go.Figure:
    """ Bar plot to visualize the count of `from` and `to` """
    order = {x: sorted(set(df[x])), y: sorted(set(df[y]))}
    df.loc[:, "txt"] = df[y] + ["*" if i else "" for i in df["is_multi"]]
    return px.bar(df, x=x, y="value", color=y, text="txt",
                  category_orders=order,
                  color_discrete_sequence=px.colors.qualitative.Light24)


def plotFromToPerSample(df: pd.DataFrame, x: str = "from",
                        title: str = "") -> list[go.Figure]:
    """ Bar plots to visualize the count of `from` and `to` and vice versa """
    df = df.rename(columns={"size": "value"})
    figs = [
        plotFromToBar(df, x, "to"),
        plotFromToBar(df[df["is_multi"] == False], x, "to"),
        # reverse
        plotFromToBar(df, "to", x),
        plotFromToBar(df[df["is_multi"] == False], "to", x),
    ]
    for fig in figs:
        fig.update_layout(title=title)

    return figs


def plotFromToPerSampleWithLevel(df: pd.DataFrame, title: str = "") -> list[go.Figure]:
    """ Bar plots including gene-level and allele-level """
    # gene-level
    df_gene = df.groupby(["from", "to", "is_multi", "correct"],
                         as_index=False).size().reset_index()
    # allele-level
    df_allele = df.groupby(["from_allele", "to", "is_multi", "correct"],
                           as_index=False).size().reset_index()
    return [
        *plotFromToPerSample(df_gene, "from", title=title),
        *plotFromToPerSample(df_allele, "from_allele", title=title),
    ]


def calcFromToStat(df: pd.DataFrame,
                   allow_multi: bool = True, method: str = "") -> pd.DataFrame:
    """ Aggregate the ratio of correct gene """
    if not allow_multi:
        df = df[df["is_multi"] == False]
    dfg = df.groupby(["to", "correct"]).size()
    dfg.name = "acc"
    dfg = dfg / dfg.groupby("to").sum()
    new_df = dfg.reset_index()
    new_df = new_df[new_df["correct"] == True]
    new_df["method"] = method
    new_df["multi"] = allow_multi
    return new_df


def plotFromToStat(df: pd.DataFrame) -> list[go.Figure]:
    """ Plot mapping stat of gene from all the samples """
    df["multi"] = df["multi"].apply(lambda i: "multi" if i else "nomulti")
    group, color = reColor(df["method"], df["multi"])
    order = {
        "multi": ["multi", "nomulti"],
        "to": sorted(set(df["to"])),
    }
    return [
        px.box(df, x="to", y="acc",
               color=group,
               color_discrete_map=color,
               category_orders=order)
          .update_layout(
              yaxis_title="Gene-level correct reads / all reads on the gene",
              title="Specificity",
          ),
        px.box(df, x="method", y="acc", color="multi",
               category_orders=order)
          .update_layout(
              yaxis_title="Specificity over all genes",
              title="Specificity",
          ),
    ]


def plotMappingFromTo() -> list[go.Figure]:
    """ The status of where the read mapped on and generated from """
    # prefix = "data6/linnil1_syn_s2022.{}.30x_s1031"
    prefix = "data6/linnil1_syn_s44.{}.30x_s444"
    cohort = [
        {"method": "hisat",           "name": f"{prefix}.index5_kir_2100_withexon_split.leftalign.mut01.graph.variant.noerrcorr.json"},
        {"method": "hisat_ab",        "name": f"{prefix}.index5_kir_2100_withexon_ab.leftalign.mut01.graph.variant.noerrcorr.json"},
        {"method": "hisat_ab_2dl1s1", "name": f"{prefix}.index5_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr.json"},
    ]
    figs = []
    df_acc = []
    for dat in cohort:
        for filename in NamePath(dat['name']).get_input_names()[:3]:
            df_from_to = calcFromToPerSample(filename)
            figs.extend(plotFromToPerSampleWithLevel(df_from_to, title=filename))
            df_acc.append(calcFromToStat(df_from_to, allow_multi=True,  method=dat["method"]))
            df_acc.append(calcFromToStat(df_from_to, allow_multi=False, method=dat["method"]))

    df = pd.concat(df_acc)
    figs.extend(plotFromToStat(df))
    return figs


def plotGeneDepth() -> list[go.Figure]:
    cohort = [{
        "method": "all",
        "name":   "data_real/hprc.{}"
                  ".index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim.variant.noerrcorr.no_multi.depth.tsv",
    }, {
        "method": "bwa_loose_filter",
        "name":   "data_real/hprc.{}.index_hs37d5.bwa.part_merge.annot_read"
                  ".index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim.variant.noerrcorr.no_multi.depth.tsv",
    }, {
        "method": "bwa_strict_filter",
        "name":   "data_real/hprc.{}.index_hs37d5.bwa.part_strict"
                  ".index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim.variant.noerrcorr.no_multi.depth.tsv",
    }]
    # for info in cohort:
    #     info['name'] = info['name'].replace("hprc.", "twbb.")

    data = []
    for info in cohort:
        for name in NamePath(info["name"]).get_input_names()[:3]:
            data.append(readSamtoolsDepth(name))
            data[-1]["method"] = info["method"]
            data[-1]["id"] = name.template_args[0]

    df = pd.concat(data)
    df["gene"] = df["gene"].str.replace("*BACKBONE", "", regex=False)

    # customze here
    df_select = df[df["gene"].isin(["KIR3DL2", "KIR3DL3"])]
    df_select = df_select.groupby(["gene", "pos", "method"], as_index=False)["depth"].median()
    print(df_select)

    group, color = reColor(df_select["method"], df_select["gene"])
    fig = px.scatter(df_select, x="pos", y="depth", color=group,
                     color_discrete_map=color, log_y=True)
    fig.update_traces(marker_size=4)
    fig.update_layout(
        legend_title_text="method - gene",
        legend_itemsizing='constant',
        xaxis_title="Positions",
        yaxis_title="Median Depth across samples",
        xaxis_type="linear",
    )
    return [fig]


if __name__ == "__main__":
    """
    extractTargetReadFromBam("data/linnil1_syn_wide.test10.00.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.no_multi.sam",
                             mapped_ref="KIR2DL1S1*BACKBONE",
                             wanted_gene="KIR2DL2")
    extractTargetReadFromBam("data/linnil1_syn_wide.test10.00.index_kir_2100_raw.mut01.hisatgenotype.errcorr.no_multi.sam",
                             mapped_ref="KIR2DL1*BACKBONE",
                             wanted_gene="KIR2DS1")
    """
    figs = []
    # figs.extend(plotBamStat())
    # figs.extend(plotGenewiseMapping())
    # figs.extend(plotMappingFromTo())
    figs.extend(plotDepth())
    print("Plot")
    app = Dash(__name__)
    app.layout = html.Div([dcc.Graph(figure=f) for f in figs])
    app.run_server(debug=True, port=8052)
