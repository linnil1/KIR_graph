# https://github.com/droeatumn/kpi
import os
from pathlib import Path
import pandas as pd
from namepipe import nt, NameTask
from kg_utils import runDocker, runShell


@nt
def setup(input_name, kpi_folder):
    if Path(kpi_folder).exists():
        return kpi_folder
    runShell(f"git clone https://github.com/droeatumn/kpi.git {kpi_folder}")
    # runShell("docker build . -t kirkpi ", cwd=kpi_folder)
    return kpi_folder


@nt
def linkSamples(input_name, data_folder):
    # 1 -> 1
    name = input_name.split('/')[0]
    output_name = os.path.join(data_folder, name + ".{}")
    new_name = output_name.format(input_name.template_args[0])
    if Path(new_name + ".read.1.fq").exists():
        return output_name
    # using hard link because we change data folder
    runShell(f"ln {input_name}.1.fq {new_name}.read.1.fq")
    runShell(f"ln {input_name}.2.fq {new_name}.read.2.fq")
    return output_name


@nt
def runKPI(input_name, kpi_folder):
    mapping_file = input_name.replace_wildcard("_merge_for_mapping")
    if Path(mapping_file + ".txt").exists():
        return input_name + ".kpi_prediction"

    with open(mapping_file + ".txt", "w") as f:
        for name in input_name.get_input_names():
            basename = Path(name).name
            print(basename + ".kpi", name + ".read.1.fq",sep="\t", file=f)
            print(basename + ".kpi", name + ".read.2.fq",sep="\t", file=f)
    folder = Path(input_name).parents[0]

    runDocker("kpi", f"./main.nf --map {mapping_file}.txt --output {folder}",
              opts=f"-v $PWD/{folder}:/opt/kpi/{folder}", workdir="/opt/kpi")
    return input_name + ".kpi_prediction"


@nt
def collectResult(input_name, kpi_folder):
    haps = pd.read_csv(f"{kpi_folder}/input/haps.txt", sep="\t")

    output_name = input_name.replace_wildcard("_merge_cn")

    cn = []
    for name in input_name.get_input_names():
        df = pd.read_csv(f"{name}.txt", sep="\t")
        haplo = df['haplotypes'][0]
        print(name, haplo)
        haplo = haplo.split("|")[0]
        a = haps[haps["nomenclature"].isin(haplo.split("+"))]
        a = a.drop(columns=["haplotype", "nomenclature", "Jiang 2012 freq", "structure"])
        a.columns = map(lambda i: "KIR" + i, a.columns)
        gene_cn = dict(a.sum(axis=0))
        gene_cn['name'] = name
        cn.append(gene_cn)

    cn = pd.DataFrame(cn)
    cn = cn.set_index("name").astype(int)
    print(cn)
    cn.to_csv(output_name + '.csv')
    return output_name


def testThreshold():
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from dash import Dash, dcc, html
    from sklearn.neighbors import KernelDensity
    from scipy.signal import argrelextrema
    import numpy as np

    # data from linnil1_syn_30x_seed87 cohert
    ori_data = np.array([57.0, 24.0, 27.0, 60.0, 29.0, 30.0, 27.0, 0.0, 27.0, 26.0, 28.0, 59.0, 60.0, 56.0, 28.0, 59.0, 0.0, 55.0, 60.0, 60.0, 0.0, 57.0, 58.0, 59.0, 60.0, 57.0, 0.0, 57.0, 49.0, 0.0, 60.0, 61.0, 30.0, 55.0, 27.0, 28.0, 26.0, 29.0, 60.0, 60.0, 56.5, 27.0, 87.0, 23.0, 28.0, 60.0, 61.0, 60.0, 27.0, 59.0, 27.0, 29.0, 59.0, 60.0, 56.0, 26.0, 28.0, 49.0, 0.0, 29.0, 29.0, 0.0, 55.0, 0.0, 28.0, 25.0, 29.0, 60.0, 60.0, 26.0, 0.0, 57.0, 49.0, 0.0, 30.0, 59.0, 31.0, 56.0, 29.0, 28.0, 24.0, 28.0, 60.0, 60.0, 26.0, 0.0, 86.0, 23.0, 27.0, 88.0, 60.0, 60.0, 28.0, 26.0, 28.0, 26.0, 29.0, 59.0, 59.0, 85.0, 55.0, 56.0, 24.0, 27.0, 59.0, 29.0, 30.0, 26.0, 0.0, 27.0, 26.0, 28.0, 59.0, 59.0, 56.0, 28.0, 145.0, 28.0, 21.0, 88.0, 122.0, 90.0, 27.0, 59.0, 0.0, 53.0, 0.0, 59.0, 60.0, 89.0, 85.0, 57.0, 49.0, 0.0, 59.0, 60.0, 30.0, 55.0, 56.0, 29.0, 0.0, 28.0, 59.0, 58.0, 48.0, 28.0])

    norm_value = 60  # 3DL3 depth
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(go.Histogram(x=ori_data / norm_value , name="Relative Depth", nbinsx=100, histnorm="probability"), secondary_y=True)
    fig.update_layout(xaxis_title="Normalized read count", yaxis_title="Normalized KDE scores")

    # GATKIR paramerts
    data = np.array(ori_data)[:, None] / norm_value
    x = np.linspace(data.min() - 0.05, data.max() + 0.05)
    kde  = KernelDensity(kernel='gaussian', bandwidth=0.075).fit(data)
    y  = kde.score_samples(x[:, None])
    mi, ma   = argrelextrema(y, np.less)[0], argrelextrema(y, np.greater)[0]
    fig.add_trace(go.Scatter(x=x, y=y / np.abs(y.min()),  name=f"GATKIR", line_color="green"))
    for cn, m in enumerate(mi):
        fig.add_vline(x=x[m], line_width=2, line_dash="dash", line_color="green", annotation_text=f"CN={cn}", annotation_font_color="green")

    # Our method
    from kg_typping_linnil1 import KDEcut
    kde = KDEcut()
    kde.fit(ori_data)
    x = np.linspace(0, 1.1, kde.points)
    y = kde.kde.score_samples(x[:, None])
    fig.add_trace(go.Scatter(x=x * kde.max / norm_value, y=y / np.abs(y.min()),  name=f"our_method", line_color="blue"))
    for cn, m in enumerate(kde.local_min):
        fig.add_vline(x=m * kde.max / norm_value, line_width=2, line_dash="dash", line_color="blue", annotation_text=f"CN={cn}", annotation_yshift=-20, annotation_font_color="blue")

    # Bad paramerts: bandwidth
    data = np.array(ori_data)[:, None] / norm_value
    x = np.linspace(data.min() - 0.05, data.max() + 0.05)
    kde  = KernelDensity(kernel='gaussian', bandwidth=0.025).fit(data)
    y  = kde.score_samples(x[:, None])
    mi, ma   = argrelextrema(y, np.less)[0], argrelextrema(y, np.greater)[0]
    fig.add_trace(go.Scatter(x=x, y=y / np.abs(y.min()),  name=f"Bad parameters", line_color="red"))
    for cn, m in enumerate(mi):
        fig.add_vline(x=x[m], line_width=2, line_dash="dash", line_color="red", annotation_text=f"CN={cn}", annotation_font_color="red")

    app = Dash(__name__)
    app.layout = html.Div(dcc.Graph(figure=fig))
    app.run_server(debug=True, port=8051)


if __name__ == "__main__":
    # testThreshold()
    kpi_folder = "kirkpi"
    kpi = None >> setup.set_args(kpi_folder)
    answer = "linnil1_syn_30x"
    answer = "linnil1_syn_30x_seed87"
    data_folder = "data3"
    Path(data_folder).mkdir(exist_ok=True)
    samples = answer + "/" + answer + ".{}.read" >> linkSamples.set_args(data_folder)
    samples >> runKPI.set_depended(-1).set_args(kpi_folder) >> collectResult.set_depended(-1).set_args(kpi_folder)
    # ignore pseduo gene in kg_eval.py
