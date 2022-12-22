import re
import subprocess
from pathlib import Path
from functools import partial
from collections import defaultdict
import numpy as np
import pandas as pd

from namepipe import compose, NameTask, BaseTaskExecutor
from graphkir.utils import runShell, getThreads
from kg_utils import runDocker, linkSamples, getAnswerFile, compareResult, buildDocker, SlurmTaskExecutor
from kg_eval import readAnswerAllele, saveCohortAllele


def extractID(name: str) -> str:
    if "hprc" in name:  # bad habit but works
        return re.findall(r"hprc\.(\w+)\.", name)[0]
    else:
        return re.findall(r"\.(\d+)\.", name)[0]


def readPingLocusCount(locus_csv: str) -> pd.DataFrame:
    ping_data = pd.read_csv(locus_csv)
    ping_data = ping_data.rename(columns={'Unnamed: 0': "sample"})
    ping_data['method'] = "PING"
    ping_data['id'] = list(map(extractID, ping_data['sample']))
    ping_data = ping_data.drop(columns=['sample'])
    return ping_data


def readAns(ans_csv):
    kir = readAnswerAllele(ans_csv)
    df = []
    for id, alleles in kir.items():
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


def cutThresholdByAns(answer_csv, sample_index):
    """ Find the threshold by answer """
    # read
    df_ping = readPingLocusCount(f"{sample_index}/locusRatioFrame.csv")
    df_ans = readAns(answer_csv)
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
    id_set = set(df[df["method"] == "PING"]["id"]) & set(df[df["method"] == "ANS"]["id"])
    df = df[df["id"].isin(id_set)]

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


def plotPing(input_name, sample_name):
    folder = str(Path(input_name).parent)
    print(folder)
    from dash import Dash, dcc, html
    import plotly.express as px

    # read data
    df_list = []
    df_list.append(readPingLocusCount(f"{folder}/locusRatioFrame.csv"))
    df_list.append(readAns(getAnswerFile(sample_name)))
    df = pd.concat(df_list)
    print(df)
    df = df.melt(["id", "method"], var_name="gene")
    df = df.sort_values(["method", "value"], ascending=[False, True])  # PING value is ascending

    # threshold
    if Path(f"{folder}/manualCopyThresholds.csv").exists():
        threshold_df = pd.read_csv(f"{folder}/manualCopyThresholds.csv")
        threshold_df.columns = ["gene", *threshold_df.columns[1:]]
        print(threshold_df)

    # plot
    figs = []
    for gene in sorted(set(df["gene"])):
        df_gene = df[df["gene"] == gene]
        fig = px.scatter(df_gene, x="id", y="value", color="method", title=gene)
        # fig = px.scatter(df_gene[df_gene['method'] == "PING"], x="id", y="value", color="method", title=gene)
        fig.update_layout(yaxis_title=f"{gene}/KIR3DL3 ratio",
                          xaxis_title="Sample ID")

        # plot threshold
        if Path(f"{folder}/manualCopyThresholds.csv").exists():
            for i in range(6):
                id = f"{i}-{i+1}"
                threshold = threshold_df[threshold_df["gene"] == gene][threshold_df.columns[i + 1]]
                if not threshold.isnull().sum() and len(threshold):
                    fig.add_hline(y=float(threshold), annotation_text=f"CN{i} - {i+1}",
                                  line_dash="dash", line_color="gray")
        figs.append(fig)

    # with open(f"{folder}/plots_of_cn_threashold.html", 'w') as f:
    #     for fig in figs:
    #         f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
    # start dash
    app = Dash(__name__)
    app.layout = html.Div(children=[dcc.Graph(figure=fig) for fig in figs])
    app.run_server(debug=True, port=8051)


def buildPing(input_name, folder):
    if Path(folder).exists():
        return folder
    runShell(f"git clone https://github.com/wesleymarin/PING.git {folder}")
    buildDocker("ping", "ping.dockerfile", folder=folder)
    return folder


def pingMain(index, folder_in, folder_out):
    runDocker("ping", f"Rscript PING_run.R", opts=f""" \
        -e INDEX={index} \
        -e RAW_FASTQ_DIR=../{folder_in} \
        -e FASTQ_PATTERN=fq \
        -e THREADS={getThreads()} \
        -e RESULTS_DIR=../{folder_out} \
    """)


def pingRun(input_name, index, sample_name):
    folder_in = str(Path(input_name).parent)
    folder_out = folder_in + ".result_"
    folder_out += index.replace("/", "_")

    if not Path(folder_out + "/manualCopyThresholds.csv").exists():
        # first time will fail
        # but we'll get locusRatioFrame.csv
        try:
            pingMain(index, folder_in, folder_out)
            # -e SHORTNAME_DELIM={threads}
        except subprocess.CalledProcessError:
            pass

        # Use answer and ratio to cut thresholds
        assert Path(folder_out + "/locusRatioFrame.csv").exists()
        answer_file = getAnswerFile(sample_name)
        print(f"Set CN threshold for PING", answer_file)
        cutThresholdByAns(answer_file, folder_out)

    # This will success
    if not Path(folder_out + "/finalAlleleCalls.csv").exists():
        pingMain(index, folder_in, folder_out)
    return folder_out + "/finalAlleleCalls"


def reorganizeResult(input_name):
    output_name = input_name + ".merge"
    if Path(output_name + ".tsv").exists():
        return output_name
    data = readPingResult(f"{input_name}.csv")
    saveCohortAllele(data, f"{output_name}.tsv")
    return output_name


def readPingResult(csv_file: str) -> dict[str, list[str]]:
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
        alleles = [i for i in alleles if "failed" not in i]
        # alleles = [i for i in alleles if "unresolved" not in i]
        called_alleles[id] = alleles
        # print(id, alleles)
    # print(called_alleles)
    return called_alleles


def removeSample(input_name, remove_name: set[str]):
    if input_name.template_args[-1] in remove_name:
        runShell(f"rm {input_name}*")
    return input_name



if __name__ == "__main__":
    samples = "linnil1_syn/linnil1_syn_s44.{}.30x_s444"
    samples = "linnil1_syn/linnil1_syn_s2022.{}.30x_s1031"
    ping_folder = f"data6/ping_{str(samples).replace('/', '_').replace('.{}', '_cohort')}"
    samples = "linnil1_syn/linnil1_syn_fakeintron1_s1214.{}.30x_s2022"
    ping_folder = f"data5/ping_{str(samples).replace('/', '_').replace('.{}', '_cohort')}"
    samples_ans = samples
    exe = BaseTaskExecutor()
    remove_sample_list = set()

    HPRC = True      # turn on HPRC samples
    TAIWANIA = False # run in TAIWANIA-HPC
    if HPRC == True:
        # requrie run kg_real.py first
        samples = "data_tmp/hprc.{}.index_hs37d5.bwa.part_strict"
        ping_folder = f"data_tmp/ping_{str(samples).replace('/', '_').replace('.{}', '_cohort')}"
        # sample_possible_ans = "data_real/hprc_merge.index_hs37d5.bwa.part_strict"
        #                       ".index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim"
        #                       ".variant.noerrcorr.no_multi.depth.p75.CNgroup_assume3DL3"
        #                       ".pv_exonfirst_1.2.compare_sum.top600"  # from data_real.py
        samples_ans = "hprc_summary"  # from kg_from_kelvin.py
        # remove_sample_list = {"HG00733", "NA19240", "HG02109", "NA21309"}
        remove_sample_list = {"HG02109", "NA21309"}  # fail by not enough KIR3DL3
        if TAIWANIA:
            # docker save localhost/linnil1/ping -o image.ping.tar.gz
            # singularity build  localhost/linnil1/ping docker-archive://image.ping.tar.gz
            exe = SlurmTaskExecutor(threads_per_sample=14, template_file="taiwania.48.template")
        else:
            samples = "data_real/ping_data_tmp_hprc_cohort.index_hs37d5.bwa.part_strict/{}"


    # main
    ping_index = compose([
        None,
        # partial(buildPing, folder="PING"),
        partial(buildPing, folder="PING20220527"),
    ])
    if not(HPRC and not TAIWANIA):
        samples = compose([
            samples,
            partial(linkSamples, data_folder=ping_folder),
            partial(removeSample, remove_name=remove_sample_list),
        ])
    ping_predict = compose([
        samples,
        NameTask(partial(pingRun, index=str(ping_index), sample_name=samples_ans), depended_pos=[-1], executor=exe),
        # partial(plotPing, sample_name=samples_ans),  # debug used: plot the gene ratio, useful when tunning the threshold
        reorganizeResult,
        partial(compareResult, sample_name=samples_ans, plot=False),
        # partial(plotPing,      sample_name=samples_ans),
    ])
    print(ping_predict)
