import asyncio
import subprocess
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor, as_completed
import re
from matplotlib.ticker import MaxNLocator
import numpy as np
from collections import defaultdict, Counter
from Bio import SeqIO
import json
import pandas as pd
from pyHLAMSA import Genemsa


def samtools(cmd, name):
    cmd = f"docker run -it --rm --name {name.replace('/', '_')} -v $PWD:/app -w /app samtools samtools {cmd} {name}"
    print(name, cmd)
    a = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    return a.stdout.decode().split("\n")


def getPairNum(name):
    """
    cmd = f"docker run -it --rm --name {name.replace('/', '_')} -v $PWD:/app -w /app samtools samtools flagstat {name}"
    print(name, cmd)
    a = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    # print(a.stdout)
    for i in a.stdout.decode().split("\n"):
    """
    num_pair = 0
    num_total = 0
    num_second = 0
    for i in samtools("flagstat", name):
        # 122050 + 0 properly paired (100.00% : N/A)
        # 122050 + 0 in total (QC-passed reads + QC-failed reads)
        # 0 + 0 secondary
        if "in total" in i:
            num_total = int(i.split()[0])
        if "secondary" in i:
            num_second = int(i.split()[0])
        if "properly paired" in i:
            num_pair = int(i.split()[0])
    return num_total - num_second, num_second, num_pair


def calc_bam_mapping():
    stat = []
    n = 10

    exc = []
    with ThreadPoolExecutor(max_workers=50) as executor:
        # tot, sec, pair = getPairNum(name)
        for i in range(n):
            name = f"data/linnil1_syn_full.{i:02d}.merge.sam"
            exc.append(executor.submit(getPairNum, name))
            # exc.append(getPairNum(name))

        for i in range(n):
            name = f"data/linnil1_syn_full.{i:02d}.split.sam"
            exc.append(executor.submit(getPairNum, name)) # exc.append(getPairNum(name))

        for i in range(n):
            name = f"data/linnil1_syn_full.{i:02d}.linear.sam"
            exc.append(executor.submit(getPairNum, name))
            # exc.append(getPairNum(name))

        for i in range(n):
            name = f"PING/PING/linnil1_syn_full_result/gc_bam_files/linnil1_syn_full.{i:02d}.read..bam"
            exc.append(executor.submit(getPairNum, name))
            # exc.append(getPairNum(name))

        for job in exc:
            # tot, sec, pair = getPairNum(name)
            tot, sec, pair = job.result()
            stat.append((tot, sec, pair))
    json.dump({'stat': stat}, open(f'tmp_bam_count.json', 'w'))


def plot_bam_mapping_fig(df, target):
    import plotly.graph_objects as go
    # Figure: Proper Paired Perecetage
    fig = go.Figure()
    methods = ['hisat-merge', 'hisat-split', 'linear', 'PING(full)']

    for i in range(len(methods)):
        fig.add_trace(go.Box(y=list(df[df['method'] == methods[i]][target]), name=methods[i]))

    for i in range(len(methods)):
        value = np.mean(df[df['method'] == methods[i]][target])
        # fig.add_annotation(text=f"{value:.3f}", x=i, y=1.05, xref="paper", yref="paper", showarrow=False)
        fig.add_annotation(text=f"mean={value:.3f}",
                           x=i, y=max(df[df['method'] == methods[i]][target]),
                           yanchor='bottom', font_size=12, showarrow=False)

    fig.update_layout(
        title="Proper Paired Percetage",
        yaxis_title="Primary paired number / Total reads (%)",
        width=800,
        height=800,
        legend=dict(
            x=0.99,
            y=0.01,
            xanchor="right",
            yanchor="bottom",
        ),
        # yaxis_range=[0.8,1],
    )
    return fig


def plotOnDash(figs):
    import dash
    from dash import dcc, html
    app = dash.Dash(__name__)
    app.layout = html.Div(children=[dcc.Graph(figure=fig) for fig in figs])
    app.run_server(port=8051, debug=True)

def plotOnDashTable(figs):
    # multiple line
    import dash
    from dash import dcc, html
    import dash_bootstrap_components as dbc
    from jupyter_dash import JupyterDash
    # app = dash.Dash(__name__)
    app = JupyterDash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
    app.layout = html.Div([
        dbc.Row([
            dbc.Col(dcc.Graph(figure=f)) for f in fig
        ])
        for fig in figs
    ])
    app.run_server(port=8051, debug=True)


def plot_bam_mapping():
    # plt.subplot(1, 3, 1)
    data = json.load(open('tmp_bam_count.json'))['stat']
    perc_pair = [pair / tot for (tot, sec, pair) in data]
    perc_secd = [sec / tot for (tot, sec, pair) in data]
    methods = ['hisat-merge', 'hisat-split', 'linear', 'PING(full)']
    n = 10

    df = pd.DataFrame(perc_pair, columns=['precision'])
    df['secondary'] = perc_secd
    df['method'] = [method for method in methods for i in range(n)]
    print(df)


    fig = plot_bam_mapping_fig(df, 'precision')
    fig.update_layout(
        title="Proper Paired Percetage",
        yaxis_title="Primary paired number / Total reads (%)",
    )

    # Figure: Secondary Paired Perecetage
    fig1 = plot_bam_mapping_fig(df, 'secondary')
    fig1.update_layout(
        title="Secondary Percetage",
        yaxis_title="Secondary amount / Total reads",
    )

    plotOnDash([fig, fig1])
    # fig.show()

    """
    import plotly.express as px
    fig = px.box(df, x="method", y="precision")
    fig.update_layout(title="Proper Paired Percetage",
                      yaxis_title="Primary paired number / Total reads (%)")

    fig = px.box(df, x="method", y="secondary")
    fig.update_layout(title="Secondary Percetage",
                      yaxis_title="Secondary amount / Total reads")
    fig.show()
    """

    """
    plt.title("Proper Paired percetage")
    plt.boxplot([perc_pair[0:n], perc_pair[n:n*2], perc_pair[n*2:n*3], perc_pair[n*3:n*4]])
    plt.ylabel("Primary paired number / Total reads (%)")
    plt.gca().set_xticklabels(methods)
    plt.savefig("tmp_bam_count_paired.png")
    plt.show()


    # plt.subplot(1, 3, 2)
    plt.title("Secondary ratio")
    plt.ylabel("Secondary amount / Total reads")
    plt.boxplot([perc_secd[n*0:n*1], perc_secd[n*1:n*2], perc_secd[n*2:n*3], perc_secd[n*3:n*4]])
    plt.gca().set_xticklabels(methods)
    plt.savefig("tmp_bam_count_secondary.png")
    plt.show()
    """

    print("Mean of perper paired percetage")
    for i in range(len(methods)):
        print(f"{methods[i]:15s} {np.mean(perc_pair[n*i:n*i+n]):.2f}")

    print("Mean of secondary ratio")
    for i in range(len(methods)):
        print(f"{methods[i]:15s} {np.mean(perc_secd[n*i:n*i+n]):.2f}")


def getMissRead(name, all_reads):
    # read all bam/sam
    miss_reads = set()
    reads = set()
    for line in samtools("view", name):
        line = line.strip()
        if not line:
            continue
        fields = line.split()
        if len(fields) < 4:
            print(line)
            continue
        reads.add(fields[0])

        # mate unmapped for read unmapped
        flag = int(fields[1])
        if (flag & 4) or (flag & 8):
            miss_reads.add(fields[0])
    # maybe the bam file did not save unmapped reads
    miss_reads.update(all_reads - reads)
    miss_count = Counter([i.split("*")[0] for i in miss_reads])
    all_read_count = Counter([i.split("*")[0] for i in all_reads])
    norm_miss_count = {gene: count / all_read_count[gene] for gene, count in miss_count.items()}
    return norm_miss_count


def missRead():
    n = 10
    sample_all_reads = []
    for i in range(n):
        # read all read id
        name = f"data/linnil1_syn_full.{i:02d}.read1.fastq"
        all_reads = set()
        for line in open(name):
            if line.startswith("@"):
                all_reads.add(line.split()[0][1:-2])
        sample_all_reads.append(all_reads)
    sample_all_reads_gene_count = [Counter([i.split("*")[0] for i in all_reads]) for all_reads in sample_all_reads]

    exc = []
    with ThreadPoolExecutor(max_workers=50) as executor:
        names = []
        names.extend([f"data/linnil1_syn_full.{i:02d}.merge.sam" for i in range(n)])
        names.extend([f"data/linnil1_syn_full.{i:02d}.split.sam" for i in range(n)])
        names.extend([f"data/linnil1_syn_full.{i:02d}.linear.sam" for i in range(n)])
        names.extend([f"PING/PING/linnil1_syn_full_result/gc_bam_files/linnil1_syn_full.{i:02d}.read..bam" for i in range(n)])
        for i, name in enumerate(names):
            # miss_reads = getMissRead(name)
            i = i % n
            exc.append(executor.submit(getMissRead, name, sample_all_reads[i]))

        sample_miss = []
        for i in exc:
            norm_miss_count = i.result()
            print(norm_miss_count)
            sample_miss.append(norm_miss_count)

    json.dump({'miss': sample_miss, 'sum': sample_all_reads_gene_count}, open(f'tmp_bam_miss.json', 'w'))


def missPlot():
    import pandas as pd
    import plotly.express as px
    n = 10
    sample_miss = json.load(open(f'tmp_bam_miss.json'))['miss']
    a = []
    methods = ['hisat-merge', 'hisat-split', 'linear', 'PING(full)']
    for i, miss in enumerate(sample_miss):
        for k, v in miss.items():
            a.append({'sample': i, 'gene': k, 'percetage': v, "method": methods[i // n]})
    df = pd.DataFrame(a)
    print(sorted(set(df['gene'])))
    fig = px.box(df, x="gene", y="percetage", color="method", category_orders={'gene': sorted(set(df['gene']))})
    fig.update_layout(title="Missed reads belong to",
                      yaxis_title="Missed_reads / total_read_counts in that gene")
    plotOnDash([fig])
    # fig.write_image("tmp_bam_miss_plot.png")
    # fig.show()


def getSecRead(name, all_reads):
    # read all bam/sam
    miss_reads = set()
    reads = set()
    for line in samtools("view", name):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) < 4:
            print(name)
            print(line)
            continue
        reads.add(fields[0])

        # mate unmapped for read unmapped
        flag = int(fields[1])
        if (flag & 256):
            miss_reads.add(fields[0])
    # maybe the bam file did not save unmapped reads
    miss_count = Counter([i.split("*")[0] for i in miss_reads])
    all_read_count = Counter([i.split("*")[0] for i in all_reads])
    norm_miss_count = {gene: count / all_read_count[gene] for gene, count in miss_count.items()}
    return norm_miss_count


def secdRead():
    n = 10
    sample_all_reads = []
    for i in range(n):
        # read all read id
        name = f"data/linnil1_syn_full.{i:02d}.read1.fastq"
        all_reads = set()
        for line in open(name):
            if line.startswith("@"):
                all_reads.add(line.split()[0][1:-2])
        sample_all_reads.append(all_reads)
    sample_all_reads_gene_count = [Counter([i.split("*")[0] for i in all_reads]) for all_reads in sample_all_reads]

    exc = []
    with ThreadPoolExecutor(max_workers=50) as executor:
        names = []
        names.extend([f"data/linnil1_syn_full.{i:02d}.merge.sam" for i in range(n)])
        names.extend([f"data/linnil1_syn_full.{i:02d}.split.sam" for i in range(n)])
        names.extend([f"data/linnil1_syn_full.{i:02d}.linear.sam" for i in range(n)])
        names.extend([f"PING/PING/linnil1_syn_full_result/gc_bam_files/linnil1_syn_full.{i:02d}.read..bam" for i in range(n)])
        for i, name in enumerate(names):
            # miss_reads = getMissRead(name)
            i = i % n
            exc.append(executor.submit(getSecRead, name, sample_all_reads[i]))

        sample_miss = []
        for i in exc:
            norm_miss_count = i.result()
            print(norm_miss_count)
            sample_miss.append(norm_miss_count)

    json.dump({'secd': sample_miss, 'sum': sample_all_reads_gene_count}, open(f'tmp_bam_secd.json', 'w'))


def secdPlot():
    import pandas as pd
    import plotly.express as px
    n = 10
    sample_miss = json.load(open(f'tmp_bam_secd.json'))['secd']
    a = []
    methods = ['hisat-merge', 'hisat-split', 'linear', 'PING(full)']
    for i, miss in enumerate(sample_miss):
        for k, v in miss.items():
            a.append({'sample': i, 'gene': k, 'percetage': v, "method": methods[i // n]})
    df = pd.DataFrame(a)
    fig = px.box(df, x="gene", y="percetage", color="method", category_orders={'gene': sorted(set(df['gene']))})
    fig.update_layout(title="Secondary reads belong to",
                      yaxis_title="Secondary_reads / total_read_counts in that gene")
    plotOnDash([fig])
    # fig.write_image("tmp_bam_secd_plot.png")
    # fig.show()


def getPrimaryRead(name):
    # read all bam/sam
    read_mapped_on = defaultdict(list)
    read_mapped_on_secd = defaultdict(list)
    for line in samtools("view", name):
    # for line in open(name):
        if line.startswith("@"):
            continue
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) < 4:
            print(name)
            print(line)
            continue
        flag = int(fields[1])
        if (flag & 2) and not (flag & 4) and not (flag & 8):
            if not (flag & 256):  # pair
                read_mapped_on[fields[0]].append(fields[2])
            else:
                read_mapped_on_secd[fields[0]].append(fields[2])

    for k, v in read_mapped_on.items():
        if len(v) != 2:
            print(k, v)
        elif v[0] != v[1]:
            print(k, v)
    for k, v in read_mapped_on_secd.items():
        if len(v) % 2:
            print(k, v)
    return read_mapped_on, read_mapped_on_secd


def primaryAcc():
    n = 10
    names = []

    # Make usre bam files didn't remove singleton and missing
    all_read_primary = []
    all_read_secondary = []
    names.extend([f"data/linnil1_syn_full.{i:02d}.split.sam" for i in range(n)])
    names.extend([f"data/linnil1_syn_full.{i:02d}.linear.sam" for i in range(n)])
    names.extend([f"PING/PING/linnil1_syn_full_result/gc_bam_files/linnil1_syn_full.{i:02d}.read..bam" for i in range(n)])
    exc = []
    with ThreadPoolExecutor(max_workers=50) as executor:
        for name in names:
            exc.append(executor.submit(getPrimaryRead, name))
            all_read_primary.append(read_primary)
            all_read_secondary.append(read_secondary)

        for i in exc:
            # read_primary, read_secondary = getPrimaryRead(name)
            read_primary, read_secondary = i.result()
            all_read_primary.append(read_primary)
            all_read_secondary.append(read_secondary)

    json.dump({'prim': all_read_primary, 'secd': all_read_secondary}, open(f'tmp_bam_mapped_on.json', 'w'))


def readAccPlot():
    print("Reading")
    import pandas as pd
    import plotly.express as px
    data = json.load(open(f'tmp_bam_mapped_on.json'))
    methods = ['hisat-split', 'linear', 'PING(full)']
    n = 10

    stat = []
    i = 0
    for read_primary, read_secondary in zip(data['prim'], data['secd']):
        print("Primary", len(read_primary))
        print("Secondary", len(read_secondary))
        tot_prim, tot_secd = 0, 0
        acc_prim, acc_secd, tot = 0, 0, 0
        for k, v in read_primary.items():
            tot_prim += len(v)
            if k.split("*")[0] == v[0].split("*")[0]:
                acc_prim += 1
        for k, v in read_secondary.items():
            tot += 1
            tot_secd += len(v)
            if (k.split("*")[0] == read_primary[k][0].split("*")[0] \
                or k.split("*")[0] in [i.split("*")[0] for i in v]):
                acc_secd += 1
        for k in read_primary.keys() - read_secondary.keys():
            tot += 1
            if k.split("*")[0] == read_primary[k][0].split("*")[0]:
                acc_secd += 1

        stat.append({
            'read': 'primary',
            'value': acc_prim / len(read_primary),
            'method': methods[i // n],
            'tot': tot_prim,
        })
        stat.append({
            'read': 'primary+secondary',
            'value': acc_secd / tot,
            'method': methods[i // n],
            'tot': tot_prim + tot_secd,
        })
        i += 1

    stat = pd.DataFrame(stat)
    print(stat.groupby(["method", "read"]).mean())

    # plot
    fig = px.box(stat, x="method", y="value", color="read")
    fig.update_layout(title="Gene-level Accuracy (per paired read)",
                      yaxis_title="Correct_reads_mapping_on / total_read_counts")

    # TODO
    """
    ('PING(full)', 'primary') 0.8766862714629227 54669.6 123
    ('PING(full)', 'primary+secondary') 0.9996833664520878 54669.6 123
    ('hisat-split', 'primary') 0.8566461611724234 54564.6 123
    ('hisat-split', 'primary+secondary') 0.9951793052557456 54564.6 123
    ('linear', 'primary') 0.8340483198512473 47049.0 123
    ('linear', 'primary+secondary') 0.9621157981826187 47049.0 123
    """
    ind_map = {'PING(full)': 1, 'hisat-split': -1, 'linear': 0, 'primary': -1, 'primary+secondary': 1}
    for i in sorted(stat.groupby(["method", "read"]).mean().itertuples()):
        fig.add_annotation(# text=f"total reads={i.tot:.0f}",
                           text=f"reads={i.tot:.0f}",
                           xref='x domain', yref='y domain',
                           x=ind_map[i.Index[0]] * 0.35 + ind_map[i.Index[1]] * 0.05 + 0.5, y=1,
                           xanchor='center', yanchor='bottom',
                           font_size=12, showarrow=False)
    plotOnDash([fig])
    # fig.write_image("tmp_bam_mapped_on.png")
    # fig.show()


# TODO: clean below codes

def plot_full_acc():
    arr_acc = []
    arr_acc_gene = []
    for i in range(10):
        tot, acc, acc_gene = 0, 0, 0
        for i in open(f"data/linnil1_syn.{i:02d}.full.sam"):
            if i and "@" == i[0]:
                continue
            row = i.split()
            if len(row) < 2:
                print(i)
            ans = row[0].split('-')[0]
            pred = row[2]
            tot += 1
            if ans == pred:
                acc += 1
            if ans.split("*")[0] == pred.split("*")[0]:
                acc_gene += 1
        arr_acc.append(acc / tot)
        arr_acc_gene.append(acc_gene / tot)

    # plt.subplot(1, 3, 3)
    plt.title("Accuracy for mapped on all kir star alleles")
    plt.boxplot([arr_acc, arr_acc_gene])
    plt.gca().set_xticklabels(['Allele', 'Gene'])
    plt.show()


def plot_each_iter():
    start_save_count = False
    loss = []
    tpm  = []
    name = []
    rec = []
    acc = []
    for i in open("tmp_out1"):
        if "Rank" in i:
            if rec:
                if len(rec) > 28:
                    # rec, tpm, name = list(zip(*sorted(zip(rec, tpm, name))[::-1]))
                    plt.suptitle(f"allele num = {len(rec)}")
                    plt.subplot(1,2,1)
                    plt.title(f"perpotion")
                    plt.bar(range(len(rec)), rec)
                    ax = plt.gca()
                    ax.set_xticks(range(len(rec)))
                    ax.set_xticklabels(name, rotation=90)
                    plt.subplot(1,2,2)
                    plt.title(f"tpm")
                    plt.bar(range(len(rec)), tpm)
                    ax = plt.gca()
                    ax.set_xticks(range(len(rec)))
                    ax.set_xticklabels(name, rotation=90)
                    plt.show()
            start_save_count = False
            rec = []
            tpm = []
            name = []
        if i.strip() and start_save_count and "Allele" in i:
            name.append(i.split()[1])
            rec.append(float(i.split()[-1]))
            tpm.append(float(i.split()[-3]))
        if i.strip() and start_save_count and "ACC" in i:
            acc.append(float(i.split()[1]))
        if "Rank 0" in i:
            loss.append(float(i.split()[-1]))
            start_save_count = True


    plt.grid()
    x = np.arange(len(loss)) + 1
    plt.xticks(x)
    plt.plot(x, loss, 'b.-')
    plt.ylabel("loss", color='b')
    ax2 = plt.gca().twinx()
    ax2.plot(x, acc, 'g.-')
    ax2.set_ylabel("Acc", color='g')
    plt.show()



def read_full_pair(fname, strict=False):
    read_map = defaultdict(list)
    for line in open(fname):
        if line and "@" == line[0]:
            continue
        if strict and "150M" not in line:
            continue
        row = line.split()
        ans = row[0].split('-')[0]
        if int(row[1]) & 64:
            read_map[row[0] + "/1"].append(row[2])
        else:
            read_map[row[0] + "/2"].append(row[2])

    pair_map = {}
    for i in set(read_map.keys()):
        if i[-1] == "1":
            id = i[:-2]
            pair_map[id] = set(read_map[id + "/1"]) & set(read_map[id + "/2"])
    return pair_map


def plotSplitAcc():
    i = 0
    pair_map = read_full_pair(f"data/linnil1_syn.{i:02d}.full.sam")
    pair_num = [len(i) for i in pair_map.values()]
    plt.title("Number of multiple alignment")
    plt.boxplot(pair_num)
    plt.show()


def plotAnsDepth():
    # name = "linnil1_syn"
    name = "linnil1_syn_full"
    # seqs = SeqIO.parse("kir_merge_sequences.fa", "fasta")
    seqs = SeqIO.parse("kir_split_full_sequences.fa", "fasta")
    seq_len = {seq.id: len(seq.seq) for seq in seqs}

    for i in range(10):
        allele_count = defaultdict(int)
        for line in open(f"{name}/{name}.{i:02d}.read..sam"):
            if line[0] == "@":
                continue
            row = line.split()
            allele_count[row[0].split("-")[0]] += 1

        alleles = [(a, b, b / seq_len[a] * 150) for a, b in allele_count.items()]
        plotDepth(sorted(alleles, key=lambda j: -j[2]))


def plotAnsGeneDepth():
    # name = "linnil1_syn"
    name = "linnil1_syn_full"
    # seqs = SeqIO.parse("kir_merge_sequences.fa", "fasta")
    seqs = SeqIO.parse("kir_split_full_backbone.fa", "fasta")
    seq_len = {seq.id: len(seq.seq) for seq in seqs}

    for i in range(10):
        allele_count = defaultdict(int)
        for line in open(f"{name}/{name}.{i:02d}.read..sam"):
            if line[0] == "@":
                continue
            row = line.split()
            allele_count[row[0].split("*")[0]] += 1

        alleles = [(a, b, b / seq_len[a + "*BACKBONE"] * 150) for a, b in allele_count.items()]
        plotDepth(sorted(alleles, key=lambda j: -j[2]))


def plotAnsDepthWithMulti():
    seqs = SeqIO.parse("kir_merge_sequences.fa", "fasta")
    seq_len = {seq.id: len(seq.seq) for seq in seqs}

    alleles_set = []
    i = 0
    name = "linnil1_syn"
    for line in open(f"{name}/{name}.{i:02d}.read..sam"):
        if line[:3] == "@SQ":
            alleles_set.append(line.split("\t")[1][3:])
    alleles_set = set(alleles_set)

    allele_count = defaultdict(int)
    pair_map = read_full_pair(f"data/{name}.{i:02d}.full.sam", strict=True)
    for id, alleles in pair_map.items():
        ans_alleles = alleles_set & alleles
        for i in ans_alleles:
            allele_count[i] += 1 / len(ans_alleles)

    alleles = [(a, b, b / seq_len[a] * 150) for a, b in allele_count.items()]
    plotDepth(sorted(alleles, key=lambda j: -j[2]))


def plotDepth(alleles, need_sort=False):
    # alleles = [(name, perpotion, tpm),]
    if need_sort:
        alleles = sorted(alleles, key=lambda i: -i[2])
    n = list(range(len(alleles)))
    plt.suptitle(f"Allele num = {len(n)}")
    allele_names, allele_count, allele_tpm = zip(*alleles)
    plt.subplot(1,2,1)
    plt.title(f"Depth")
    plt.bar(n, allele_tpm)
    ax = plt.gca()
    ax.set_xticks(n)
    ax.set_xticklabels(allele_names, rotation=90)
    plt.subplot(1,2,2)
    plt.title(f"perpotion")
    plt.bar(n, allele_count)
    ax = plt.gca()
    ax.set_xticks(n)
    ax.set_xticklabels(allele_names, rotation=90)
    plt.show()


def evaluateHisatMapPlot():
    from pipeline_hisat_kir import HisatTyping
    i = 0
    names = []
    for i in range(1, 10):
        names.append(f"data/linnil1_syn.0{i}.merge.pair.tmp")
    # names = names[:1]
    typ = HisatTyping()

    acc_and_arr, acc_sum_arr = [], []

    for name in names:
        # typ.mainPerSample(name)
        typ.name = name
        typ.gene = "KIR"
        typ.readBam()
        acc_and, acc_sum = typ.evaluateHisatMap()
        acc_and_arr.append(acc_and)
        acc_sum_arr.append(acc_sum)

    plt.title("Allele Accuracy for Hisat2 mapping")
    plt.boxplot([acc_and_arr, acc_sum_arr])
    plt.ylabel("Accuracy")
    plt.gca().set_xticklabels(['AND', 'SUM'])
    plt.show()


def checkSecondIsPair():
    i = 0
    def getName():
        for i in range(10):
            yield f"PING/PING/linnil1_syn_full_result/gc_bam_files/linnil1_syn_full.{i:02d}.read..bam"
            yield f"data/linnil1_syn_full.{i:02d}.split.sam"
            yield f"data/linnil1_syn_full.{i:02d}.linear.sam"

    for name in getName():
        # read all bam/sam
        secd_read = 0
        reads = 0
        bad_secd = 0
        reads_id = set()
        for line in samtools("view", name):
            line = line.strip()
            if line.startswith("@"):
                continue
            if not line:
                continue
            fields = line.split()
            if len(fields) < 4:
                print(line)
                continue
            flag = int(fields[1])
            reads += 1
            reads_id.add(fields[0])

            if (flag & 4) or (flag & 8):
                continue
            if (flag & 256):
                # four flag is proper pair with 256(secdonary)
                if flag not in [355, 403, 339, 419]:
                    # There is still a very little 337 353
                    # print(line)
                    bad_secd += 1
                secd_read += 1
        print(len(reads_id), reads, secd_read, bad_secd)


def mappingRatePerStage():
    def getBam():
        for i in range(1):
            yield {
                'method': "ans",
                'sample': f"linnil1_syn_full/linnil1_syn_full.{i:02d}.read..sam",
            }
            """
            yield {
                'method': "hisat-noab",
                'sample': f"data/linnil1_syn_full.{i:02d}.noab.pair.bam",
            }
            yield {
                'method': "hisat-noab-typing",
                'sample': f"data/linnil1_syn_full.{i:02d}.noab.pair.tmp.sam",
            }
            """
            yield {
                'method': "hisat-noab-sec",
                'sample': f"data/linnil1_syn_full.{i:02d}.noab.sec_pair.bam",
            }
            yield {
                'method': "hisat-noab-typing",
                'sample': f"data/linnil1_syn_full.{i:02d}.noab.sec_pair.tmp.sam",
            }
            yield {
                'method': "hisat-noab-typing-sec-noNH",
                'sample': f"data/linnil1_syn_full.{i:02d}.noab.sec_pair.noNH.tmp.sam",
            }
            yield {
                'method': "hisat-noab-typing-sec-noNH-avg",
                'sample': f"data/linnil1_syn_full.{i:02d}.noab.sec_pair.noNH.tmp.sam",
                'weighted': True,
            }
            """
            yield {
                'method': "hisat-merge",
                'sample': f"data/linnil1_syn_full.{i:02d}.merge.pair.bam",
            }
            yield {
                'method': "hisat-merge-typing",
                'sample': f"data/linnil1_syn_full.{i:02d}.merge.pair.tmp.sam"
            }
            yield {
                'method': "hisat-merge-typing-noNH",
                'sample': f"data/linnil1_syn_full.{i:02d}.merge.pair.noNH.tmp.group.sam"
            }
            yield {
                'method': "hisat-split",
                'sample': f"data/linnil1_syn_full.{i:02d}.split.pair.bam",
            }
            yield {
                'method': "hisat-split-typing",
                'sample': f"data/linnil1_syn_full.{i:02d}.split.pair.tmp.sam",
            }
            yield {
                'method': "hisat-split-typing-noNH",
                'sample': f"data/linnil1_syn_full.{i:02d}.split.pair.noNH.tmp.sam",
            }
            """

    def countGene(f):
        return Counter(map(lambda i: i.split('*')[0], filter(lambda i: i.strip(), samtools("view", f))))

    def countGeneWithWeight(f):
        count_gene = defaultdict(int)
        lines = filter(lambda i: i.strip(), samtools("view", f))
        for i in lines:
            NH = re.findall("NH\:i\:(\d+)", i)
            if NH:
                weight = 1/int(NH[0])
            else:
                weight = 1
            count_gene[i.split('*')[0]] += weight # shoudl /2 bcz pair but fine
        return count_gene

    num_per_file = []
    for data in getBam():
        if data.get("weighted"):
            data.update(**countGeneWithWeight(data['sample'].replace("full", "wide")))
        else:
            data.update(**countGene(data['sample'].replace("full", "wide")))
        # data.update(**countGene(data['sample']))
        num_per_file.append(data)

    df = pd.DataFrame(num_per_file)
    df = df.drop(columns=['weighted'])
    df = pd.melt(df, id_vars=["method", "sample"], 
                 var_name="gene", value_name="read_count")
    print(df)


    df1 = df
    df['percetage'] = df.apply(lambda i: float(i['read_count']) / float(df[ (df['method'] == "ans")  & (df['gene'] == i['gene'])]['read_count']), axis=1)

    # plot
    import plotly.express as px
    # color should be same for each catelog
    # colors = px.colors.qualitative.Dark24
    # color_map = { method: colors[i] for i, method in enumerate(sorted(set(df["method"]))) }

    figs = []
    figs.append(px.bar(df, x="gene", y="read_count", color="method", barmode="group", category_orders={'gene': sorted(set(df['gene']))}))
    # df1 = df[df['method'] != "ans"]
    figs.append(px.box(df, x="gene", y="percetage", color="method", category_orders={'gene': sorted(set(df['gene']))}))
    plotOnDash(figs)


def diffBetweenAllele(msa, title):
    import plotly.express as px
    figs = []
    bs = msa.get_variantion_base()
    print(title)
    print("Total length", msa.get_length())
    print("Total base diff", len(bs))
    print("Continuous base diff", len(re.findall(r"gDNA\s+(\d+)", msa.format_alignment_diff())))
    df = pd.DataFrame(bs, columns = ['pos'])
    figs.append( px.histogram(df, x='pos', title=title, width=600, height=800) )
    figs.append( px.histogram(df, x='pos', title=title, width=600, height=800, histnorm='probability').update_layout(yaxis_tickformat = '.2%') )
    return figs


def plotDiffBetweenPair():
    figs = []
    msa = Genemsa.load_msa(f"kir_merge_full.save.fa", f"kir_merge_full.save.gff")

    # gene = "KIR2DS1"
    # base_gene = "KIR2DL1"
    # gene = "KIR2DL2"
    # base_gene = "KIR2DL3"
    gene = "KIR2DS3"
    base_gene = "KIR2DS5"

    # All genes
    # for g in set(map(lambda i: i.split("*")[0], msa.alleles.keys())):
    for g in [gene, base_gene]:
        submsa = msa.select_allele(f"({g})|({base_gene})").shrink()
        figs.append(diffBetweenAllele(submsa, f"{base_gene} vs {g}"))

    # diff variant of two = all variant - variant1 - variant2
    submsa = msa.select_allele(f"({gene})|({base_gene})").shrink()
    submsa_base0 = submsa.get_variantion_base()
    submsa_base1 = submsa.select_allele(f"({gene})").get_variantion_base()
    submsa_base2 = submsa.select_allele(f"({base_gene})").get_variantion_base()
    bases = set(submsa_base0) - set(submsa_base1) - set(submsa_base2)

    # plot diff base
    """
    merged_bases = []
    right = 5
    for b in sorted(bases):
        if len(merged_bases) and merged_bases[-1][1] + right * 2 >= b:
            merged_bases[-1][1] = b
        else:
            merged_bases.append([b, b])
    for b_left, b_right in merged_bases:
        print(submsa.format_alignment_from_center(b_left, right=b_right - b_left + right))
    """

    # plot
    title = f"Different between {base_gene} and {gene} (exclude internal variant)"
    print(title, len(bases))
    import plotly.express as px
    # .update_layout(yaxis_range=[0,70])
    df = pd.DataFrame(bases, columns = ['pos'])
    figs.append([px.histogram(df, x='pos', title=title, width=600, height=800),
                 px.histogram(df, x='pos', title=title, width=600, height=800, histnorm='probability').update_layout(yaxis_tickformat = '.2%') ])
    plotOnDashTable(figs)

def getBase(msa, reg):
    return set(msa.select_allele(reg).get_variantion_base())

def plotDiffBetweenGene():
    figs = []
    msa = Genemsa.load_msa(f"kir_merge_full.save.fa", f"kir_merge_full.save.gff")
    del msa.alleles['KIR*BACKBONE']
    del msa.alleles['KIR2DS4*0010103']
 
    # remove 2DL5 short sequences
    for i in map(lambda i: i[0], filter(lambda i: len(i[1].replace("-", "")) < 3000, msa.select_allele("KIR2DL5A.*").alleles.items())):
        print(f"delete {i}")
        del msa.alleles[i]
    for i in map(lambda i: i[0], filter(lambda i: len(i[1].replace("-", "")) < 3000, msa.select_allele("KIR2DL5B.*").alleles.items())):
        print(f"delete {i}")
        del msa.alleles[i]
 

    genes = sorted(set(map(lambda i: i.split("*")[0], msa.alleles.keys())))
    # genes.remove("KIR3DL3")
    genes.remove("KIR3DP1")
    # genes.remove("KIR2DL5A")
    # genes.remove("KIR2DL5B")
    gene_id = {gene: id for id, gene in enumerate(genes)}

    bases = {}
    for gene in genes:
        bases[gene] = set(msa.select_allele(f"{gene}").get_variantion_base())

    diff_genes = []

    exes = {}
    with ProcessPoolExecutor(max_workers=4) as executor:
        # run concurrent
        for gene1 in genes:
            for gene2 in genes:
                if gene1 != gene2:
                    exes[(gene1, gene2)] = executor.submit(getBase, msa, f"({gene1})|({gene2})")

        # gene1 vs gene2
        for gene1 in genes:
            for gene2 in genes:
                if gene1 != gene2:
                    print(gene1, gene2)
                    b = exes[(gene1, gene2)].result()
                    diff_genes.append({'from': gene1, 'to': gene2, 'value': len(b - bases[gene1] - bases[gene2])})
                else:
                    diff_genes.append({'from': gene1, 'to': gene2, 'value': len(bases[gene1])})
    diff_genes_df = pd.DataFrame(diff_genes)

    # pandas to 2d array
    data = np.zeros((len(genes), len(genes)))
    for d in diff_genes:
        data[gene_id[d['from']], gene_id[d['to']]] = d['value']

    import plotly.express as px
    import plotly.graph_objects as go
    figs.append(px.imshow(data, text_auto=True, color_continuous_scale='RdBu_r',
                          width=1200, height=1200,
                          labels=dict(color="Base"),
                          x=list(gene_id.keys()), y=list(gene_id.keys())
                ).update_xaxes(side="top")
    )

    # dash cannot show text on image
    # figs[0].write_image("genes_diff_all.png")
    # figs[0].write_image("genes_diff.png")  # no 3DP1
    figs[0].write_image("genes_diff_noshort.png")  # no 3DP1 + no 2DL5 short
    plotOnDash(figs)

    # exc.append(executor.submit(self.batchInsertVariant, data, check_exist))
    # for res in tqdm(as_completed(exc), total=len(exc)):
    #     insert_count, total = res.result()

if __name__ == "__main__":
    # calc_bam_mapping()
    # plot_bam_mapping()
    # missRead()
    # missPlot()
    # secdRead()
    # secdPlot()
    # primaryAcc()
    # readAccPlot()
    # checkSecondIsPair()
    # mappingRatePerStage()
    # plotDiffBetweenGene()
    # plotDiffBetweenPair()

    # tmp
    # evaluateHisatMapPlot()
    # plot_full_acc()
    # plotAnsDepth() 
    # plotAnsGeneDepth() 
    # plotAnsDepthWithMulti() 
    pass
