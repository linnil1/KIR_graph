import re
import subprocess
from pathlib import Path
from functools import partial
from collections import Counter, defaultdict

from dash import Dash, dcc, html
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from namepipe import (
    NameTask,
    compose,
    BaseTaskExecutor,
    ConcurrentTaskExecutor,
    NamePath,
)
from graphkir.utils import (
    getThreads,
    getGeneName,
    setThreads,
)
from kg_utils import runShell, runDocker
from graphkir.external_tools import setEngine
from graphkir.plot import showPlot
from kir.kir_pipe import FileMod, Executor
from kir.sakauekir import SakaueKir
from kir.t1k import T1k
from kir.kpi import KPI
from kir.ping import PING

from kg_utils import linkSamples, getAnswerFile, compareResult, SlurmTaskExecutor
from kg_eval import readAnswerAllele


class NamePipeFileMod(FileMod):
    """Using namepipe as our pipeline filename handler"""

    def getID(self, name: NamePath) -> str:
        return name.template_args[0]

    def listFiles(self, name: NamePath) -> list[NamePath]:
        return name.list_names()

    def replaceWildcard(self, name: NamePath, new_name: str) -> NamePath:
        return name.replace_wildcard(new_name)


class MainExecutor(Executor):
    """Connect the Executor from graphkir/utils"""

    def runDocker(
        self,
        image: str,
        cmd: str,
        cwd: str | None = None,
        opts: str = "",
    ) -> subprocess.CompletedProcess[str]:
        """run docker container"""
        return runDocker(image, cmd, cwd=cwd, opts=opts)


def createSakaueKirPloidy(input_name: NamePath, sample_name: NamePath) -> NamePath:
    """Read answer allele as input SakaueKir's CN result"""
    output_name = input_name + ".answer_cn"
    if Path(output_name + ".csv").exists():
        return output_name

    ans = readAnswerAllele(getAnswerFile(sample_name))

    def renameGene(i: str) -> str:
        return {
            "KIR2DS3": "KIR2DS3;KIR2DS5",
            "KIR2DS5": "KIR2DS3;KIR2DS5",
            "KIR2DL5A": "KIR2DL5A;KIR2DL5B",
            "KIR2DL5B": "KIR2DL5A;KIR2DL5B",
        }.get(i, i)

    ids = []
    dfs = []
    for id, alleles in ans.items():
        ids.append(id)
        dfs.append(Counter(map(renameGene, map(getGeneName, alleles))))

    df = pd.DataFrame(dfs)
    df = df.set_axis(ids, axis=0)
    df = df.T.fillna(0).astype(int)
    df.to_csv(output_name + ".csv")
    print(df)
    return output_name


def sakauekirRun(
    samples: NamePath, data_folder: str, use_answer: bool = True
) -> NameTask:
    """Run SakaueKir"""
    gatkir = SakaueKir(
        threads=getThreads(),
        file_adapter=NamePipeFileMod,
        executor=MainExecutor,
    )
    folder = compose(
        [
            None,
            gatkir.download,
        ]
    )
    index = folder.output_name

    samples_bam = compose(
        [
            samples,
            partial(linkSamples, data_folder=data_folder),
            partial(gatkir.bwa, index=index),
            gatkir.addGroup,
            gatkir.markDupliate,
        ]
    )
    samples_bam >> partial(gatkir.deepVariant, index=index)

    if use_answer:
        samples_cn = compose(
            [
                samples_bam,
                partial(createSakaueKirPloidy, sample_name=samples),
            ]
        )
    else:
        samples_cn = compose(
            [
                samples_bam,
                partial(gatkir.analysisTK, index=index),
                gatkir.getCoverage,
                NameTask(gatkir.ploidyEstimate, depended_pos=[-1]),
            ]
        )

    samples_call = compose(
        [
            samples_bam,
            partial(gatkir.beforeHC, ploidy_name=samples_cn.output_name),
            partial(gatkir.haplotypeCaller, index=index),
            NameTask(
                gatkir.jointGenotype, depended_pos=[-1], func_kwargs=dict(index=folder)
            ),
            gatkir.beforeCalling,
            partial(gatkir.calling, index=index),
            NameTask(gatkir.mergeCalling, depended_pos=[-1]),
        ]
    )
    compose(
        [
            samples_call,
            NameTask(gatkir.mergeResult, depended_pos=[-1]),
            partial(compareResult, sample_name=samples),
        ]
    )
    compose(
        [
            samples_call,
            NameTask(partial(gatkir.mergeResult, select_all=True), depended_pos=[-1]),
            partial(compareResult, sample_name=samples),
        ]
    )
    return samples_call


def t1kRun(samples: NamePath, data_folder: str) -> NameTask:
    """Run T1K"""
    t1k = T1k(
        threads=getThreads(),
        file_adapter=NamePipeFileMod,
        executor=MainExecutor,
    )
    index = compose(
        [
            None,
            t1k.download,
            partial(t1k.build, ipd_version="2100"),
        ]
    )
    samples = compose(
        [
            samples,
            partial(linkSamples, data_folder=data_folder),
            partial(t1k.run, index=index.output_name),
            NameTask(t1k.mergeResult, depended_pos=[-1]),
            partial(compareResult, sample_name=samples),
        ]
    )
    return samples


def kpiRun(samples: NamePath, data_folder: str) -> NameTask:
    """Run kpi"""
    kpi = KPI(
        threads=getThreads(),
        file_adapter=NamePipeFileMod,
        executor=MainExecutor,
    )
    index = kpi.download()
    samples = compose(
        [
            samples,
            partial(linkSamples, data_folder=data_folder),
            NameTask(kpi.run, func_kwargs=dict(index=index), depended_pos=[-1]),
            NameTask(kpi.mergeResult, func_kwargs=dict(index=index), depended_pos=[-1]),
            partial(compareResult, sample_name=samples),
        ]
    )
    return samples


def pingRun(
    samples: NamePath,
    data_folder: str,
    version: str = "20220527",
    remove_sample_list: set[str] = set(),
    use_slurm: bool = False,
    answer_name: str = "",
    show_plot_and_break: bool = False,
) -> NameTask:
    """Run ping"""
    ping = PING(
        threads=getThreads(),
        version=version,
        file_adapter=NamePipeFileMod,
        executor=MainExecutor,
    )
    if not answer_name:
        answer_name = samples
    index = ping.download()

    samples = compose(
        [
            samples,
            partial(linkSamples, data_folder=data_folder),
            NameTask(ping.migrateSample, depended_pos=[-1]),
            partial(pingRemoveSample, remove_name=remove_sample_list),
        ]
    )
    folder_out = ping.getOutputFolder(samples.output_name, index)

    if show_plot_and_break:
        showPlot(pingPlot(folder_out, answer_name))
        exit()

    # change executor
    if use_slurm:
        executor = SlurmTaskExecutor(
            threads_per_sample=14, template_file="research/taiwania.14.template"
        )
        ping.setThreads(14)
    else:
        executor = BaseTaskExecutor()

    # If you break at this line. Fine
    # Try again after editing manualCopyNumberFrame.csv
    try:
        samples = samples >> NameTask(
            ping.main, func_kwargs=dict(index=index), executor=executor
        )
    except subprocess.CalledProcessError:
        pingPredictCNByAnswer(folder_out, answer_name)
        samples = samples >> NameTask(
            ping.main, func_kwargs=dict(index=index), executor=executor
        )

    ping.setThreads(getThreads())
    samples >> partial(ping.mergeResult, use_novel=True)
    samples = compose(
        [
            samples,
            partial(ping.mergeResult, use_novel=False),
            partial(compareResult, sample_name=answer_name),
        ]
    )
    return samples


def pingPredictCNByAnswer(folder_out: str, sample_name: str, save: bool = True):
    """Predict the gene-depth-ratio threshold"""
    df_ping = PING.readGeneDepthRatio(f"{folder_out}/locusRatioFrame.csv")
    df_ans = readAnswerGeneCN(getAnswerFile(sample_name))
    """
    id method  KIR2DL1  KIR2DL2  KIR2DL4  KIR2DL5
    0  00    ans      0.5      1.0      1.0      1.5
    1  01    ans      1.0      0.5      1.0      1.0
    2  02    ans      1.0      1.0      2.0      1.0
    3  03    ans      0.5      1.0      1.5      0.5
    """
    df = pd.concat([df_ping, df_ans], ignore_index=True)
    df = df.melt(["id", "method"], var_name="gene")
    df = df.sort_values(
        ["method", "value"], ascending=[False, True]
    )  # PING value is ascending
    """
    id method     gene     value
    0    00   PING  KIR3DP1  0.343889
    1    01   PING  KIR3DP1  0.344004
    """
    id_set = set(df[df["method"] == "PING"]["id"]) & set(
        df[df["method"] == "ANS"]["id"]
    )
    df = df[df["id"].isin(id_set)]
    print("ignore " + str(set(df["id"]) - id_set))

    # cut threshold and plot
    threshold_list = []
    for gene in sorted(set(df["gene"])):
        df_part = df[df["gene"] == gene]
        threshold = pingCalcThreshold(
            df_part[df["method"] == "ANS"], df_part[df["method"] == "PING"]
        )
        threshold_dict = {"gene": gene}
        for i, c in enumerate(threshold):
            threshold_dict[f"{i}-{i+1}"] = c
        threshold_list.append(threshold_dict)

    columns = ["gene"] + [f"{i}-{i+1}" for i in range(6)]
    df_threshold = pd.DataFrame(threshold_list)
    df_threshold = df_threshold[df_threshold["gene"] != "KIR3DL3"]
    df_threshold = df_threshold.reindex(columns=columns)
    df_threshold = df_threshold.fillna("NA")
    if save:
        df_threshold.to_csv(f"{folder_out}/manualCopyThresholds.csv", index=False)
    return df, df_threshold


def readAnswerGeneCN(ans_csv: str) -> pd.DataFrame:
    """Calculate copy number from answer allele"""
    kir = readAnswerAllele(ans_csv)
    df_list = []
    for id, alleles in kir.items():
        d: dict[str, float] = defaultdict(float)
        for allele in alleles:
            gene_name = allele[:7]  # 2DL5A 2DL5B are in the same group
            d[gene_name] += 1 / 2  # CN of 3DL3 = 2
        df_list.append(
            {
                **d,
                "id": id,
                "method": "ANS",
            }
        )
    df = pd.DataFrame(df_list).fillna(0)
    return df


def pingRemoveSample(folder: str, remove_name: set[str]) -> NamePath:
    for name in Path(folder).iterdir():
        if Path(name).name.split(".")[1] in remove_name:  # pattern: id.{id}.read.xx
            runShell(f"rm {name}*")
    return folder


def pingPlot(folder: str, sample_name: str = "") -> list[go.Figure]:
    """Plot ping's CN graph"""
    # read ping depth data
    df_list = []
    df_ping = PING.readGeneDepthRatio(f"{folder}/locusRatioFrame.csv")
    df_list.append(df_ping)

    # read answer
    if sample_name:
        df_list.append(readAnswerGeneCN(getAnswerFile(sample_name)))

    # merge ping depth and real depth
    df = pd.concat(df_list)
    print(df)
    df = df.melt(["id", "method"], var_name="gene")
    df = df.sort_values(
        ["method", "value"], ascending=[False, True]
    )  # PING value is ascending

    # read threshold
    # if Path(f"{folder}/manualCopyThresholds.csv").exists():
    df_threshold = pd.read_csv(f"{folder}/manualCopyThresholds.csv")
    df_threshold = df_threshold.set_axis(["gene", *df_threshold.columns[1:]], axis=1)
    print(df_threshold)

    # plot
    figs = []
    for gene in sorted(set(df["gene"])):
        df_gene = df[df["gene"] == gene]
        fig = px.scatter(df_gene, x="id", y="value", color="method", title=gene)
        # fig = px.scatter(df_gene[df_gene['method'] == "PING"], x="id", y="value", color="method", title=gene)
        fig.update_layout(yaxis_title=f"{gene}/KIR3DL3 ratio", xaxis_title="Sample ID")

        # plot threshold
        for i in range(6):
            id = f"{i}-{i+1}"
            threshold = df_threshold[df_threshold["gene"] == gene][
                df_threshold.columns[i + 1]
            ]
            if not threshold.isnull().sum() and len(threshold):
                fig.add_hline(
                    y=float(threshold),
                    annotation_text=f"CN{i} - {i+1}",
                    line_dash="dash",
                    line_color="gray",
                )
        figs.append(fig)

    # with open(f"{folder}/plots_of_cn_threashold.html", 'w') as f:
    #     for fig in figs:
    #         f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
    # start dash
    return figs


def pingCalcThreshold(df_ans: pd.DataFrame, df_ping: pd.DataFrame) -> list[float]:
    """
    Cut the threshold by comparing answer and real gene depth ratio

    Answer: 0    0.5  0.5   1.5  1.5
    count:  0.1  0.2  0.21  0.4  0.5

    Return:   0.15      0.3 0.3    0.6
    CN_explain:  0   1     2    3      4
    """
    ans_count = list(np.round(df_ans["value"] * 2).astype(np.int_))
    ping_count = list(df_ping["value"])
    now_cn = 0  # for ans
    prev_count = 0  # for ping
    threshold = []

    for i in range(len(ping_count)):
        while ans_count[i] != now_cn:
            now_cn += 1
            threshold.append((prev_count + ping_count[i]) / 2)
        prev_count = ping_count[i]
    threshold.append(prev_count + 0.5)
    return threshold


def testKDEThreshold():
    from plotly.subplots import make_subplots
    from sklearn.neighbors import KernelDensity
    from scipy.signal import argrelextrema
    from graphkir.cn_model import KDEcut
    from graphkir.plot import showPlot

    # data from linnil1_syn_30x_seed87 cohert
    ori_data = np.array([57.0, 24.0, 27.0, 60.0, 29.0, 30.0, 27.0, 0.0, 27.0, 26.0, 28.0, 59.0, 60.0, 56.0, 28.0, 59.0, 0.0, 55.0, 60.0, 60.0, 0.0, 57.0, 58.0, 59.0, 60.0, 57.0, 0.0, 57.0, 49.0, 0.0, 60.0, 61.0, 30.0, 55.0, 27.0, 28.0, 26.0, 29.0, 60.0, 60.0, 56.5, 27.0, 87.0, 23.0, 28.0, 60.0, 61.0, 60.0, 27.0, 59.0, 27.0, 29.0, 59.0, 60.0, 56.0, 26.0, 28.0, 49.0, 0.0, 29.0, 29.0, 0.0, 55.0, 0.0, 28.0, 25.0, 29.0, 60.0, 60.0, 26.0, 0.0, 57.0, 49.0, 0.0, 30.0, 59.0, 31.0, 56.0, 29.0, 28.0, 24.0, 28.0, 60.0, 60.0, 26.0, 0.0, 86.0, 23.0, 27.0, 88.0, 60.0, 60.0, 28.0, 26.0, 28.0, 26.0, 29.0, 59.0, 59.0, 85.0, 55.0, 56.0, 24.0, 27.0, 59.0, 29.0, 30.0, 26.0, 0.0, 27.0, 26.0, 28.0, 59.0, 59.0, 56.0, 28.0, 145.0, 28.0, 21.0, 88.0, 122.0, 90.0, 27.0, 59.0, 0.0, 53.0, 0.0, 59.0, 60.0, 89.0, 85.0, 57.0, 49.0, 0.0, 59.0, 60.0, 30.0, 55.0, 56.0, 29.0, 0.0, 28.0, 59.0, 58.0, 48.0, 28.0])

    norm_value = 60  # 3DL3 depth
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(
        go.Histogram(
            x=ori_data / norm_value,
            name="Relative Depth",
            nbinsx=100,
            histnorm="probability",
        ),
        secondary_y=True,
    )
    fig.update_layout(
        xaxis_title="Normalized read count", yaxis_title="Normalized KDE scores"
    )

    # GATKIR paramerts
    data = np.array(ori_data)[:, None] / norm_value
    x = np.linspace(data.min() - 0.05, data.max() + 0.05)
    kde = KernelDensity(kernel="gaussian", bandwidth=0.075).fit(data)
    y = kde.score_samples(x[:, None])
    mi, ma = argrelextrema(y, np.less)[0], argrelextrema(y, np.greater)[0]
    fig.add_trace(
        go.Scatter(x=x, y=y / np.abs(y.min()), name=f"GATKIR", line_color="green")
    )
    for cn, m in enumerate(mi):
        fig.add_vline(
            x=x[m],
            line_width=2,
            line_dash="dash",
            line_color="green",
            annotation_text=f"CN={cn}",
            annotation_font_color="green",
        )

    # Our method
    kde = KDEcut()
    kde.fit(ori_data)
    x = np.linspace(0, 1.1, kde.points)
    y = kde.kde.score_samples(x[:, None])
    fig.add_trace(
        go.Scatter(
            x=x * kde.x_max / norm_value,
            y=y / np.abs(y.min()),
            name=f"our_method",
            line_color="blue",
        )
    )
    for cn, m in enumerate(kde.local_min):
        fig.add_vline(
            x=m * kde.x_max / norm_value,
            line_width=2,
            line_dash="dash",
            line_color="blue",
            annotation_text=f"CN={cn}",
            annotation_yshift=-20,
            annotation_font_color="blue",
        )

    # Bad paramerts: bandwidth
    data = np.array(ori_data)[:, None] / norm_value
    x = np.linspace(data.min() - 0.05, data.max() + 0.05)
    kde = KernelDensity(kernel="gaussian", bandwidth=0.025).fit(data)
    y = kde.score_samples(x[:, None])
    mi, ma = argrelextrema(y, np.less)[0], argrelextrema(y, np.greater)[0]
    fig.add_trace(
        go.Scatter(x=x, y=y / np.abs(y.min()), name=f"Bad parameters", line_color="red")
    )
    for cn, m in enumerate(mi):
        fig.add_vline(
            x=x[m],
            line_width=2,
            line_dash="dash",
            line_color="red",
            annotation_text=f"CN={cn}",
            annotation_font_color="red",
        )

    showPlot([fig])


def testPingWithKDE():
    from graphkir.cn_model import KDEcut

    folder = "data3/ping_linnil1_syn_30x_seed87.result"

    df_ping = PING.readGeneDepthRatio(f"{folder}/locusRatioFrame.csv")
    df_threshold = pd.read_csv(f"{folder}/manualCopyThresholds.csv")
    df_threshold = df_threshold.set_axis(["gene", *df_threshold.columns[1:]], axis=1)

    df = pd.concat([df_ping])
    print(df)
    print(df_threshold)
    df = df.melt(["id", "method"], var_name="gene")
    df = df.sort_values(
        ["method", "value"], ascending=[False, True]
    )  # PING value is ascending

    # plot
    figs = []
    for gene in sorted(set(df["gene"])):
        df_gene = df[df["gene"] == gene]
        fig = px.scatter(df_gene, x="id", y="value", title=gene)
        fig.update_layout(yaxis_title=f"{gene}/KIR3DL3 ratio", xaxis_title="Sample ID")
        figs.append(fig)
        kde = KDEcut()
        kde.bandwidth = 0.08
        kde.fit(df_gene["value"])
        figs.extend(kde.plot())
        figs[-1].update_layout(title=f"{gene}'s copy number estimation by KDE")

    showPlot(figs)


if __name__ == "__main__":
    # testKDEThreshold()
    # testPingWithKDE()
    # exit()
    data_folder = "data"
    # samples = "linnil1_syn/linnil1_syn_s44.{}.30x_s444"
    samples = "linnil1_syn/linnil1_syn_s2022.{}.30x_s1031"
    # samples = "linnil1_syn/linnil1_syn_s2022.{}.15x_s1031"
    setThreads(2)
    NameTask.set_default_executor(ConcurrentTaskExecutor(threads=10))
    use_slurm = False
    t1kRun(samples, data_folder)
    sakauekirRun(samples, data_folder, use_answer=True)
    setThreads(25)
    kpiRun(samples, data_folder)
    # exit()

    setThreads(25)
    remove_sample_list: set[str] = set()
    answer_name = ""
    # HPRC sample on 
    # samples = "data_real/hprc.{}.index_hs37d5.bwa.part_strict"
    # data_folder = "data_real"
    # answer_name = "hprc_summary"  # from kg_from_kelvin.py  and uncomment show_plot_and_break
    # remove_sample_list = {"HG02109", "NA21309"}

    TAIWANIA = False
    if TAIWANIA:  # Run this in TAIWANIA (for hprc sample)
        setEngine("singularity_linnil1")
        data_folder = "data_tmp"
        samples = data_folder + "/" + Path(samples).name
        use_slurm = True

    pingRun(
        samples,
        data_folder,
        version="wgs",  # 20220527, wgs
        remove_sample_list=remove_sample_list,
        use_slurm=use_slurm,
        answer_name=answer_name,
        # show_plot_and_break=True,
    )
