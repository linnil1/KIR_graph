from dash import Dash, dcc, html, Input, Output
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import os
import json
import subprocess
from collections import defaultdict, Counter


def plotOnDash(figs):
    app = Dash(__name__)
    app.layout = html.Div([dcc.Graph(figure=f) for f in figs])
    app.run_server(debug=True, port=8051)


def samtools(cmd, name):
    cmd = f"docker run -it --rm --name {name.replace('/', '_')} -v $PWD:/app -w /app samtools samtools {cmd} -@30 {name}"
    print(cmd)
    a = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    return a.stdout.decode().split("\n")


def getStat(filename):
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
        'total': num_total,
        'secd': num_second,
        'pair': num_pair,
    }


def plot_bam_mapping():
    filename = "bam_map_stat"

    if not os.path.exists(filename + ".csv"):
        # Figure: Proper Paired Perecetage
        data = [
            {'method': "answer",    'id': 0, 'file': "linnil1_syn_wide/linnil1_syn_wide.00.read..sam"},
            {'method': "raw_split", 'id': 0, 'file': "data/linnil1_syn_wide.00.kir_2100_raw.mut01.bam"},
            {'method': "merge",     'id': 0, 'file': "data/linnil1_syn_wide.00.kir_2100_merge.mut01.nosingle.bam"},
            {'method': "linear",    'id': 0, 'file': "data/linnil1_syn_wide.00.kir_2100_raw_cons.bowtie2.sam"},
            {'method': "full",      'id': 0, 'file': "data/linnil1_syn_wide.00.kir_2100_raw_full.bowtie2.sam"},
            {'method': "ping",      'id': 0, 'file': "data/linnil1_syn_wide.00.kir_2100_raw_full.ping.sam"},
        ]
        
        for i in data:
            i.update(getStat(i['file']))
        df = pd.DataFrame(data)
        ans = df[ (df['method'] == 'answer') ]
        ans = dict(zip(ans['id'], ans['total']))
        df['ans_total'] = df.apply(lambda i: ans[i.id], axis=1)
        df['pair_perc'] = df['pair'] / df['ans_total']
        df['secd_perc'] = df['secd'] / df['ans_total']
        # df.to_csv(filename + ".csv", index=False)
    else:
        df = pd.read_csv(filename + ".csv")

    fig0 = px.box(df, x="method",  y="pair_perc", title="Proper Paired Proportion",
                  labels={'pair_perc': "Primary paired reads / Total reads"})
    fig1 = px.box(df, x="method",  y="secd_perc", title="Secondary Proportion",
                  labels={'secd_perc': "Secondary reads / Total reads"})
    return [fig0, fig1]


def getMappedPos(filename):
    data = defaultdict(list)
    for i in samtools("view", filename):
        fields = i.split("\t")
        if len(fields) < 3:
            print(i)
            continue
        pos = -1
        if fields[3] != "*":
            pos = int(fields[3])
        # flag  backbone  position
        data[fields[0]].append(
            (int(fields[1]), fields[2].split('-')[0], pos)
        )
    return { 'mapped_pos': data }


def plot_missing():
    filename = "bam_mapping_detail"

    if not os.path.exists(filename + ".json"):
        # Figure: Proper Paired Perecetage
        data = [
            {'method': "answer",    'id': 0, 'file': "linnil1_syn_wide/linnil1_syn_wide.00.read..sam"},
            {'method': "raw_split", 'id': 0, 'file': "data/linnil1_syn_wide.00.kir_2100_raw.mut01.bam"},
            {'method': "merge",     'id': 0, 'file': "data/linnil1_syn_wide.00.kir_2100_merge.mut01.nosingle.bam"},
            {'method': "linear",    'id': 0, 'file': "data/linnil1_syn_wide.00.kir_2100_raw_cons.bowtie2.sam"},
            # {'method': "full",   'id': 0, 'file': "data/linnil1_syn_wide.00.kir_2100_raw_full.bowtie2.sam"},
        ]
        
        for i in data:
            i.update(getMappedPos(i['file']))
        # print(data)
        json.dump(data, open(filename + ".json", "w"))
    else:
        data = json.load(open(filename + ".json"))

    def findAnswer(id):
        for sample in data:
            if sample['method'] == "answer" and sample['id'] == id:
                return sample

    gene_missing_data = []
    for sample in data:
        print(sample['method'])
        if sample['method'] == "answer":
            continue
        reads0 = findAnswer(sample['id'])['mapped_pos']
        reads1 = sample['mapped_pos']

        # primary is not paired = miss
        reads1 = {read_name: mapped_pos for read_name, mapped_pos in reads1.items() \
                  if len(list(filter(lambda i: i[1] != '*' and not (i[0] & 256), mapped_pos))) >=2}
        total = Counter([name.split("*")[0] for name in reads0])
        miss  = Counter([name.split("*")[0] for name in set(reads0.keys()) - set(reads1.keys())])
        # for i in set(reads0.keys()) - set(reads1.keys()):
        #     print(reads0[i])
        print(total, miss)

        # accuracy
        primary_acc = defaultdict(int)
        for read_name, mapped_pos in sample['mapped_pos'].items():
            gene = read_name.split("*")[0]
            for flag, ref, pos in mapped_pos:
                # not secondary, not unmapped, not in proper paired
                if flag & 256 or flag & 4 or flag & 8 or not(flag & 2):
                    continue
                if gene == ref.split("*")[0]:
                    primary_acc[gene] += 1
        for gene in primary_acc:
            primary_acc[gene] = primary_acc[gene] / total[gene] / 2

        for gene in set(total.keys()):
            gene_missing_data.append({
                'id': sample['id'],
                'method': sample['method'],
                'gene': gene,
                'miss_percetage': miss.get(gene, 0) / total[gene],
                'miss_count': miss.get(gene, 0),
                'gene_accuracy_primary': primary_acc[gene],
            })


    df = pd.DataFrame(gene_missing_data)
    print(df)
    gene_order = {'gene': sorted(set(df['gene']))}
    fig0 = px.box(df, x="gene", y="miss_count", color="method", category_orders=gene_order)
    fig0.update_layout(title="Missed reads belong to",
                       yaxis_title="Missed_reads")
    fig1 = px.box(df, x="gene", y="miss_percetage", color="method", category_orders=gene_order)
    fig1.update_layout(title="Missed reads belong to (Proportion)",
                       yaxis_title="Missed_reads / total_read_counts in that gene")
    fig2 = px.box(df, x="gene", y="gene_accuracy_primary", color="method", category_orders=gene_order)
    fig2.update_layout(title="Gene-level accuracy",
                        yaxis_title="Correct read / total read (%)")
    fig3 = px.box(df, x="method", y="gene_accuracy_primary", category_orders=gene_order)
    fig3.update_layout(title="Gene-level accuracy",
                        yaxis_title="Correct read / total read %")

    return [fig0, fig1, fig2, fig3]


if __name__ == '__main__':
    figs = []
    figs.extend(plot_bam_mapping())
    figs.extend(plot_missing())
    plotOnDash(figs)

