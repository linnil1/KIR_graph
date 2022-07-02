import os
import subprocess
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd

from namepipe import nt, NameTask

from kg_eval import EvaluateKIR, extractID
from kg_utils import runDocker, runShell, threads


def readPingLocusCount(locus_csv):
    ping_data = pd.read_csv(locus_csv)
    ping_data = ping_data.rename(columns={'Unnamed: 0': "sample"})
    ping_data['method'] = "PING"
    ping_data['id'] = list(map(extractID, ping_data['sample']))
    ping_data = ping_data.drop(columns=['sample'])
    return ping_data


def readAns(ans_csv):
    kir = EvaluateKIR(ans_csv)
    df = []
    for id, alleles in kir.ans.items():
        d = defaultdict(int)
        d['id'] = id
        d['method'] = "ANS"
        for allele in alleles:
            gene_name = allele[:7]  # 2DL5A 2DL5B are in the same group
            d[gene_name] += 1 / 2  # CN of 3DL3 = 2
        df.append(d)
    df = pd.DataFrame(df).fillna(0)
    return df


def calcThreshold(df_ans, df_ping):
    """
    Answer: 0    0.5  0.5   1.5  1.5
    count:  0.1  0.2  0.21  0.4  0.5

    Return:   0.15      0.3 0.3    0.6
    CN_explain:  0   1     2    3      4
    """
    ans_count  = list(np.round(df_ans['value'] * 2).astype(np.int))
    ping_count = list(df_ping['value'])
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


def saveThreshold(threshold_list, filename):
    """ Save to PING's manualCopyThresholds.csv format """
    columns = ["gene"] + [f"{i}-{i+1}" for i in range(6)]
    df = pd.DataFrame(threshold_list)
    df = df[df["gene"] != "KIR3DL3"]
    df = df.reindex(columns=columns)
    df = df.fillna("NA")
    df.to_csv(filename, index=False)
    print(df)
    return df


def cutThresholdByAns(ans_csv, sample_index):
    """ Find the threshold by answer """
    # read
    df_ping = readPingLocusCount(f"{sample_index}/locusRatioFrame.csv")
    df_ans = readAns(ans_csv)
    """
       id method  KIR2DL1  KIR2DL2  KIR2DL4  KIR2DL5  
       0  00    ans      0.5      1.0      1.0      1.5
       1  01    ans      1.0      0.5      1.0      1.0
       2  02    ans      1.0      1.0      2.0      1.0
       3  03    ans      0.5      1.0      1.5      0.5
    """
    df = pd.concat([df_ping, df_ans], ignore_index=True)
    df = df.melt(["id", "method"], var_name="gene")
    df = df.sort_values(["method", "value"], ascending=[False, True])  # PING value is ascending
    """
         id method     gene     value
         0    00   PING  KIR3DP1  0.343889
         1    01   PING  KIR3DP1  0.344004
    """

    # cut threshold and plot
    threshold_list = []
    for gene in sorted(set(df["gene"])):
        df_part = df[df["gene"] == gene]
        threshold = calcThreshold(df_part[df["method"] == "ANS"], df_part[df["method"] == "PING"])
        threshold_dict = {'gene': gene}
        for i, c in enumerate(threshold):
            threshold_dict[f"{i}-{i+1}"] = c
        threshold_list.append(threshold_dict)

    # save and plot
    threshold_df = saveThreshold(threshold_list, f"{sample_index}/manualCopyThresholds.csv")
    return df, threshold_df


def plotPing():
    from dash import Dash, dcc, html
    import plotly.express as px

    # main
    name = "linnil1_syn_30x_seed87"
    name = "linnil1_syn_30x"
    ratio_df, threshold_df = cutThresholdByAns(
        f"{name}/{name}.summary.csv",
        f"PING/data_{name}.result"
    )

    # plot
    figs = []
    for gene in sorted(set(ratio_df["gene"])):
        df_part = ratio_df[ratio_df["gene"] == gene]
        threshold = threshold_df[threshold_df["gene"] == gene].to_dict('records')
        if not threshold:
            continue
        threshold = threshold[0]
        print(threshold)
        fig = px.scatter(df_part, x="id", y="value", color="method", title=gene)
        for i in range(6):
            id = f"{i}-{i+1}"
            if threshold.get(id) != "NA":
                fig.add_hline(y=threshold[id], annotation_text=f"CN{i} - {i+1}",
                              line_dash="dash", line_color="gray")
        figs.append(fig)

    # start dash
    app = Dash(__name__)
    app.layout = html.Div(children=[dcc.Graph(figure=fig) for fig in figs])
    app.run_server(debug=True)


@nt
def linkSamples(input_name, data_folder):
    # 1 -> 1
    name = input_name.split('/')[0]
    output_name = os.path.join(data_folder, name + ".{}")
    new_name = output_name.format(input_name.template_args[0])
    if Path(new_name + ".read.1.fq").exists():
        return output_name
    # using hard link because we change data folder
    runShell(f"ln -s ../{input_name}.1.fq {new_name}.read.1.fq")
    runShell(f"ln -s ../{input_name}.2.fq {new_name}.read.2.fq")
    return output_name


@nt
def buildPing(folder):
    if os.path.exists(folder):
        return folder
    runShell(f"git clone https://github.com/wesleymarin/PING.git {folder}")
    runShell(f"docker build {folder} -t ping")
    runShell(f"ln -s ../{folder} {folder}")
    return folder


@nt
def pingCopyFile(input_name):
    # generate isolated fastq folder
    name = Path(input_name).name
    output_folder = f"{Path(input_name).parents[0]}/ping_{name.split('.')[0]}"
    Path(output_folder).mkdir(exist_ok=True)
    if Path(f"{output_folder}/{name}.read.1.fq").exists():
        return output_folder
    runShell(f"ln -s ../{name}.read.1.fq {output_folder}/{name}.read.1.fq")
    runShell(f"ln -s ../{name}.read.2.fq {output_folder}/{name}.read.2.fq")
    return output_folder


def pingMain(index, folder_in, folder_out):
    assert index == "PING"
    runDocker("ping", f"Rscript {index}/PING_run.R", opts=f""" \
        -e RAW_FASTQ_DIR=../{folder_in} \
        -e FASTQ_PATTERN=fq \
        -e THREADS={threads} \
        -e RESULTS_DIR=../{folder_out} \
    """)


@nt
def ping(input_name, index):
    assert index == "PING"
    folder_in = input_name
    folder_out = input_name + ".result"

    # first time will fail
    # but we'll get locusRatioFrame.csv
    if not os.path.exists(folder_out + "/locusRatioFrame.csv"):
        try:
            pingMain(index, folder_in, folder_out)
            # -e SHORTNAME_DELIM={threads}
        except subprocess.CalledProcessError:
            pass

        # Use answer and ratio to cut thresholds
        assert os.path.exists(folder_out + "/locusRatioFrame.csv")
        # data2_linnil1_syn_30x_seed87/
        # -> linnil1_syn_30x_seed87
        name = input_name.split('/')[-1].split('.')[0].split("_", maxsplit=1)[1]
        print(f"Set CN threshold fro PING", name)
        cutThresholdByAns(
            f"{name}/{name}.summary.csv",
            folder_out
        )

    # This will success
    if not os.path.exists(folder_out + "/finalAlleleCalls.csv"):
        pingMain(index, folder_in, folder_out)
    return folder_out + "/finalAlleleCalls"


@nt
def pingResult(input_name, answer):
    kir = EvaluateKIR(f"{answer}/{answer}.summary.csv")
    ping_called_alleles = readPingResult(f"{input_name}.csv")
    kir.compareCohert(ping_called_alleles)
    return input_name


def readPingResult(csv_file):
    """ Read ping result """
    # read
    data = pd.read_csv(csv_file, sep=',')
    data = data.rename(columns={"Unnamed: 0": "name"})
    """
    print(data)
    name                      KIR3DP1                                        KIR2DS35
    linnil1_syn_wide.00.read. KIR3DP1*026+KIR3DP1*null                       KIR2DS3*009+KIR2DS3*024+KIR2DS5*02701
    linnil1_syn_wide.01.read. KIR3DP1*00302+KIR3DP1*03201 KIR3DP1*00304+KIR3
    """
    called_alleles = {}
    for sample in data.to_dict('records'):
        id = extractID(sample['name'])
        alleles = []
        for gene, alleles_str in sample.items():
            if gene == "name":
                continue
            alleles.extend(alleles_str.split(' ')[0].split('+'))
        alleles = [i for i in alleles if "null" not in i]
        # alleles = [i for i in alleles if "unresolved" not in i]
        called_alleles[id] = alleles
        print(id, alleles)
    return called_alleles


if __name__ == "__main__":
    answer = "linnil1_syn_30x_seed87"
    data_folder = "data3"
    Path(data_folder).mkdir(exist_ok=True)
    samples = f"{answer}/{answer}" + ".{}.read" >> linkSamples.set_args(data_folder)

    ping_index = "PING" >> buildPing
    ping_predict = samples >> pingCopyFile >> ping.set_args(index=str(ping_index)) 
    ping_predict >> pingResult.set_args(answer=answer)
    print(ping_predict)