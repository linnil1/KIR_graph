from dash import Dash, dcc, html, Input, Output
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import os
import json
import subprocess
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor


def get7Name(id):
    return id[:7]


def getGeneName(id):
    return id.split("*")[0]


def getAlleleName(id):
    return id.split("-")[0]


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
        id = 0
        data = [[
            {'method': "answer",           'id': id, 'file': f"linnil1_syn_wide/linnil1_syn_wide.{id:02d}.read..sam"},
            {'method': "linear",           'id': id, 'file': f"data/linnil1_syn_wide.{id:02d}.kir_2100_raw_cons.bowtie2.sam"},
            {'method': "full",             'id': id, 'file': f"data/linnil1_syn_wide.{id:02d}.kir_2100_raw_full.bowtie2.sam"},
            # {'method': "ping",             'id': 0, 'file': f"data/linnil1_syn_wide.00.kir_2100_raw_full.ping.sam"},
            {'method': "hisat_merge",      'id': id, 'file': f"data/linnil1_syn_wide.{id:02d}.kir_2100_merge.mut01.bam"},
            {'method': "hisat_merge_type", 'id': id, 'file': f"data/linnil1_syn_wide.{id:02d}.kir_2100_merge.mut01.hisatgenotype.sam"},
            {'method': "hisat_raw",        'id': id, 'file': f"data/linnil1_syn_wide.{id:02d}.kir_2100_raw.mut01.bam"},
            {'method': "hisat_raw_type",   'id': id, 'file': f"data/linnil1_syn_wide.{id:02d}.kir_2100_raw.mut01.hisatgenotype.sam"},
            {'method': "hisat_splitab",    'id': id, 'file': f"data/linnil1_syn_wide.{id:02d}.kir_2100_ab.mut01.bam"},
        ] for id in range(5)]
        data = [j for i in data for j in i]
        print(data)

        # for i in data:
        #     i.update(getStat(i['file']))
        with ProcessPoolExecutor() as executor:
            for d, result in zip(data, executor.map(getStat, map(lambda i: i['file'], data))):
                d.update(result)

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


def getMappedOn(filename):
    data = defaultdict(list)
    for i in samtools("view", filename):
        fields = i.split("\t")
        if len(fields) < 3:
            print(i)
            continue
        # read_id: (flag  backbone)
        data[fields[0]].append(
            (int(fields[1]), fields[2].split('-')[0])
        )
    return {
        'mapping': data,
        'total': len(data),
    }


def custom_get_missing(total, reads, **kwargs):
    data = {gene: {'total': num, 'count': 0, 'miss': 0} for gene, num in total.items()}

    for read_name, mapping_info in reads.items():
        data[getGeneName(read_name)]['count'] += 1
        if len(list(filter(lambda i: i[1] != '*' and not (i[0] & 256), mapping_info))) <2:
            data[getGeneName(read_name)]['miss'] += 1

    for gene, d in data.items():
        d['miss'] += d['total'] - d['count']  # maybe removed by mapper
        yield {
            'gene': gene,
            'total': d['total'],
            'count': d['count'],
            'miss': d['miss'],
            'miss_perc': d['miss'] / d['total'],
        }


def custom_get_proper_mapped(total, reads, **kwargs):
    data = {gene: {'total': num, 'count': 0, 'pair': 0} for gene, num in total.items()}

    for read_name, mapping_info in reads.items():
        data[getGeneName(read_name)]['count'] += 1
        if len(list(filter(lambda i: (i[0] & 2), mapping_info))) >= 2:
            data[getGeneName(read_name)]['pair'] += 1

    for gene, d in data.items():
        yield {
            'gene': gene,
            'total': d['total'],
            'count': d['count'],
            'pair_num': d['pair'],
            'pair_perc': d['pair'] / d['total'],
        }


def custom_secondary_count(total, reads, **kwargs):
    data = {gene: {'total': num, 'count': 0, 'secd': 0} for gene, num in total.items()}
    for read_name, mapping_info in reads.items():
        data[getGeneName(read_name)]['count'] += 1
        data[getGeneName(read_name)]['secd'] += len(list(filter(lambda i: (i[0] & 256), mapping_info)))

    for gene, d in data.items():
        yield {
            'gene': gene,
            'total': d['total'],
            'count': d['count'],
            'secd_num': d['secd'],
            'secd_perc': d['secd'] / d['count'],
        }


def custom_acc(total, reads, getGeneNameWrap):
    data = {gene: {'total': num, 'count': 0, 'primary': 0, 'secondary': 0} for gene, num in total.items()}

    for read_name, mapping_info in reads.items():
        data[getGeneName(read_name)]['count'] += 1
        # remove non-pair
        mapping_info = list(filter(lambda i: (i[0] & 2), mapping_info))

        primary = list(filter(lambda i: i[1] != '*' and not (i[0] & 256), mapping_info))
        if primary and getGeneNameWrap(primary[0][1]) == getGeneNameWrap(read_name):
            data[getGeneName(read_name)]['primary'] += 1

        second = list(filter(lambda i: i[1] != '*', mapping_info))
        if second and any(getGeneNameWrap(read[1]) == getGeneNameWrap(read_name) for read in second):
            data[getGeneName(read_name)]['secondary'] += 1

    for gene, d in data.items():
        yield {
            'gene': gene,
            'total': d['total'],
            'correct': d['primary'],
            'acc': d['primary'] / d['count'],
            'type': 'primary',
        }
        yield {
            'gene': gene,
            'total': d['total'],
            'correct': d['secondary'],
            'acc': d['secondary'] / d['count'],
            'type': 'secondary',
        }


def plot_bam_custom_evalute():
    filename = "results/bam_mapping_detail"
    filename += "_1"

    if not os.path.exists(filename + ".json"):
        # Figure: Proper Paired Perecetage
        id = 0
        data = []
        for id in range(10):
            data.extend([
                {'method': "answer",    'id': id, 'same_ab': False, 'file': f"linnil1_syn_wide/linnil1_syn_wide.{id:02d}.read..sam"},
                {'method': "raw_split", 'id': id, 'same_ab': True, 'file': f"data/linnil1_syn_wide.{id:02d}.kir_2100_raw.mut01.bam"},
                {'method': "ab_split",  'id': id, 'same_ab': False, 'file': f"data/linnil1_syn_wide.{id:02d}.kir_2100_ab.mut01.bam"},
                {'method': "merge",     'id': id, 'same_ab': True, 'file': f"data/linnil1_syn_wide.{id:02d}.kir_2100_merge.mut01.bam"},
                {'method': "linear",    'id': id, 'same_ab': True, 'file': f"data/linnil1_syn_wide.{id:02d}.kir_2100_raw_cons.bowtie2.sam"},
                # {'method': "full",   'id': 0, 'file': "data/linnil1_syn_wide.00.kir_2100_raw_full.bowtie2.sam"},
            ])

        # for i in data:
        #     i.update(getMappedOn(i['file']))
        with ProcessPoolExecutor() as executor:
            for d, result in zip(data, executor.map(getMappedOn, map(lambda i: i['file'], data))):
                d.update(result)

        # print(data)
        json.dump(data, open(filename + ".json", "w"))
    else:
        print('Reading')
        data = json.load(open(filename + ".json"))

    def reColor(arr1, arr2):
        group_name = arr1 + " - " + arr2
        colors = (px.colors.qualitative.Dark2, px.colors.qualitative.Set2)

        name_to_color = {}
        for i, m in enumerate(set(arr1)):
            for j, t in enumerate(set(arr2)):
                name_to_color[m + " - " + t] = colors[j][i]
        return group_name, name_to_color

    def findAnswerSample(id):
        for sample in data:
            if sample['method'] == "answer" and sample['id'] == id:
                return sample

    def runAllFiles(data, stat_func):
        custom_stats = []
        for sample in data:
            if sample['method'] == "answer":
                continue
            print(sample['method'], sample['id'])
            reads_ans = findAnswerSample(sample['id'])['mapping']
            reads_file = sample['mapping']
            gene_ans_total = Counter(map(lambda name: getGeneName(name), reads_ans.keys()))

            if sample['same_ab']:
                stats = stat_func(gene_ans_total, reads_file, getGeneNameWrap=lambda i: i[:7])
            else:
                stats = stat_func(gene_ans_total, reads_file, getGeneNameWrap=getGeneName)

            custom_stats.extend([{
                'method': sample['method'],
                'id': sample['id'],
                **stat,
            } for stat in stats])

        df = pd.DataFrame(custom_stats)
        print(df)
        return df

    figs = []

    # missing
    df = runAllFiles(data, custom_get_missing)
    gene_order = {'gene': sorted(set(df['gene']))}
    figs.append(
        px.box(df, x="gene", y="miss_perc", color="method", category_orders=gene_order)
          .update_layout(title="Missed reads belong to (Proportion)",
                         yaxis_title="Missed_reads / total_read_counts in that gene")
    )
    figs.append(
        px.box(df, x="method", y="miss_perc", color="method", category_orders=gene_order)
          .update_layout(title="Missed reads Per Method (Proportion)",
                         yaxis_title="Missed_reads / total_read_counts in that gene")
    )

    # Any proper reads
    df = runAllFiles(data, custom_get_proper_mapped)
    gene_order = {'gene': sorted(set(df['gene']))}
    figs.append(
        px.box(df, x="gene", y="pair_perc", color="method", category_orders=gene_order)
          .update_layout(title="Proper mapped",
                         yaxis_title="Correct read / total read (%)")
    )
    figs.append(
        px.box(df, x="method", y="pair_perc", category_orders=gene_order)
          .update_layout(title="Proper mapped",
                         yaxis_title="Pair-mapped read / total read (%)")
    )

    # secondary
    df = runAllFiles(data, custom_secondary_count)
    gene_order = {'gene': sorted(set(df['gene']))}
    figs.append(
        px.box(df, x="gene", y="secd_perc", color="method", category_orders=gene_order)
          .update_layout(title="Secondary count",
                         yaxis_title="secondary read / total read (%)")
    )
    figs.append(
        px.box(df, x="method", y="secd_perc", category_orders=gene_order)
          .update_layout(title="Secondary Count",
                         yaxis_title="secondary read / total read (%)")
    )

    # gene-acc
    df = runAllFiles(data, custom_acc)
    gene_order = {'gene': sorted(set(df['gene']))}
    group, color = reColor(df["method"], df["type"])
    figs.append(
        px.box(df, x="gene", y="acc", color=group,
               category_orders=gene_order,
               color_discrete_map=color,
              )
          .update_layout(title="Gene-level accuracy(Primary)",
                         yaxis_title="Correct read / total read (%)")
    )
    figs.append(
        px.box(df, x="method", y="acc", color="type", category_orders=gene_order)
          .update_layout(title="Gene-level accuracy Per Method (Proportion)",
                         yaxis_title="Correct read / total read (%)")
    )
    return figs


if __name__ == '__main__':
    figs = []
    # figs.extend(plot_bam_mapping())
    figs.extend(plot_bam_custom_evalute())
    app = Dash(__name__)
    app.layout = html.Div([dcc.Graph(figure=f) for f in figs])
    app.run_server(debug=True, port=8051)
