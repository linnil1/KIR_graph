"""
Plot utility for statistic (Not for accuracy)
"""
from typing import Iterable
from concurrent.futures import ProcessPoolExecutor

from Bio import SeqIO
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from .kir_cn import readSamtoolsDepth
from .cn_model import loadCNModel
from .utils import runDocker, runShell


def plotCN(filename_json: str) -> list[go.Figure]:
    """ Load CN mdoel and plot """
    dist = loadCNModel(filename_json)
    name = filename_json.rsplit(".", 1)[0]
    return dist.plot(title=name)


def plotGeneDepths(file_samtools_depths: str, title: str = "") -> list[go.Figure]:
    """ Plot depth of each gene """
    df = readSamtoolsDepth(file_samtools_depths)
    df["gene"] = df["gene"].str.replace("*BACKBONE", "", regex=False)
    fig_line = px.line(df, x="pos", y="depth", color="gene",
                       title=f"Read Depth {title}")
    fig_line.update_xaxes(type='linear')
    fig_box = px.box(df, y="gene", x="depth",
                     orientation="h",  # horizontal
                     title=f"Read Depth {title}")
    fig_box.update_xaxes(type='linear')
    return [fig_line, fig_box]


def readSamtoolsFlagstat(bamfile: str) -> dict[str, int]:
    """ read samtools flag stat """
    num_pair = 0
    num_total = 0
    num_second = 0
    proc = runDocker("samtools", f"samtools flagstat -@4 {bamfile}",
                     capture_output=True)
    for i in proc.stdout.split("\n"):
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
        'total': num_total,
        'secd': num_second,
        'pair': num_pair,
    }


def readFastqID(fastq_file: str) -> list[str]:
    """ Read read ID in fastq """
    seqs = SeqIO.parse(fastq_file, "fastq")
    return list(map(lambda i: str(i.id), seqs))


def plotReadMappingStat(bam_files: list[str],
                        fastq_files: list[str] = [],
                        methods: list[str] | None = None) -> list[go.Figure]:
    """
    Plot read mapping result

    Parameters:
      bam_files: the path of bamfile
      methods:
        The corresponding method for the bamfile.
        Leave None if only one method.
    """
    with ProcessPoolExecutor() as executor:
        stats = executor.map(readSamtoolsFlagstat, bam_files)
        if fastq_files:
            reads = executor.map(readFastqID, fastq_files)
    runShell("stty echo opost")

    df = pd.DataFrame(list(stats))
    df['name'] = bam_files

    if methods is None:
        df['method'] = "test"
    else:
        df['method'] = methods

    if fastq_files:
        df['read_total'] = [len(i) * 2 for i in reads]
    else:
        df['read_total'] = df['total']

    print(df)
    df['pair_perc'] = df['pair'] / df['read_total']
    df['secd_perc'] = df['secd'] / df['read_total']
    print(df)

    fig0 = px.box(df, x="method",  y="pair_perc", title="Proper Paired Ratio",
                  labels={'pair_perc': "Primary paired reads / Total reads"})
    fig1 = px.box(df, x="method",  y="secd_perc", title="Secondary Ratio",
                  labels={'secd_perc': "Secondary reads / Total reads"})
    return [fig0, fig1]


def showPlot(figs: Iterable[go.Figure]) -> None:
    """ Show all the figure in website default localhost:8051 """
    from dash import Dash, dcc, html
    app = Dash(__name__)
    app.layout = html.Div([dcc.Graph(figure=f) for f in figs])
    app.run_server(debug=True, port=8051)


def savePlot(html_filename: str, figs: Iterable[go.Figure]) -> None:
    """ Save all the figure in html_filename """
    with open(html_filename, 'w') as f:
        for fig in figs:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
