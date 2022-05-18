from collections import defaultdict
import numpy as np
import pandas as pd
from kg_eval import EvaluateKIR, extractID


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
    # df.to_csv(filename, index=False)
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


if __name__ == "__main__":
    from dash import Dash, dcc, html
    import plotly.express as px

    # main
    ratio_df, threshold_df = cutThresholdByAns(
        "linnil1_syn_30x/linnil1_syn_30x.summary.csv",
        "PING/data_linnil1_syn_30x.result"
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
