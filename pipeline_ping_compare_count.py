"""
Compare ping gene-count vs answer
"""
import pandas as pd
from collections import Counter
import plotly.express as px
import dash
from dash import dcc, html


def plotOnDash(figs):
    app = dash.Dash(__name__)
    app.layout = html.Div(children=[dcc.Graph(figure=fig) for fig in figs])
    app.run_server(debug=True)


def readCount():
    # Ping
    ping_data = pd.read_csv("PING/PING_test/linnil1_syn_wide_result/locusRatioFrame.csv")
    ping_data = ping_data.rename(columns={'Unnamed: 0': "sample"})
    ping_data['method'] = "PING"
    ping_data['sample'] = ping_data['sample'].str[:-6]

    # Answer
    real_data = pd.read_csv("PING/PING_test/linnil1_syn_wide/summary.csv", sep="\t", names=['sample', 'haplo', 'alleles'])
    real_data['method'] = "ANS"
    # string to {gene: couont}
    for i, d in real_data.iterrows():
        allele = d['alleles'].split("_")
        allele = map(lambda i: i.split("*")[0], allele)
        gene_counts = Counter(allele)
        for gene_name, gene_count in gene_counts.items():
            real_data.loc[i, gene_name] = gene_count
    real_data = real_data.fillna(0)
    real_data['KIR2DL5'] = real_data['KIR2DL5A'] + real_data['KIR2DL5B']
    real_data = real_data.drop(columns=['KIR2DL5A', 'KIR2DL5B'])
    real_data = real_data.drop(columns=['alleles','haplo'])

    all_data = pd.concat([ping_data, real_data])
    return all_data


def geneOrder(all_data, gene):
    pindex = all_data['method'] == "PING"
    rank = list(all_data[pindex][gene].rank(method='min'))
    rank_set = set()
    for i, num in enumerate(rank):
        while num in rank_set:
            num += 1
        rank_set.add(num)
        rank[i] = num
    all_data.loc[pindex, 'order'] = rank
    for _, data in all_data[pindex].iterrows():
        all_data.loc[ (all_data['sample'] == data['sample']) & (all_data['method'] == "ANS"), 'order' ] = data['order'] - 0.1
    all_data.loc[(all_data['method'] == "ANS"), gene ] /= 2
    return all_data


all_datas = readCount()
print(all_datas.columns)

figs = []
for gene in sorted(set(all_datas) - set(['sample', 'method'])):
    # gene = "KIR2DL1"
    print(gene)
    all_data = geneOrder(all_datas.copy(), gene)
    fig = px.scatter(all_data, x="order", y=gene, color="method", title=gene)
    figs.append(fig)
plotOnDash(figs)
