from dash import Dash, dcc, html, Input, Output
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
import pandas as pd
import os
import json
import subprocess
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor
from kg_utils import runDocker, samtobam


def get7Name(id):
    return id[:7]


def getGeneName(id):
    return id.split("*")[0]


def getAlleleName(id):
    return id.split("-")[0]


def samtools(cmd, name):
    proc = runDocker("samtools", f"samtools {cmd} -@4 {name}",
                     capture_output=True)
    return proc.stdout.split("\n")


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


def plotBamStat():
    # sample_index = "data/linnil1_syn_wide.test10"
    sample_index = "data/linnil1_syn_30x"
    filename = f"{sample_index}.bam_stat"
    filename = f"{sample_index}.bam_stat_vg"
    answer_name = "linnil1_syn_wide/linnil1_syn_wide"
    answer_name = "linnil1_syn_30x/linnil1_syn_30x"

    if not os.path.exists(filename + ".json"):
        # Figure: Proper Paired Perecetage
        id = 0
        data = []
        for id in range(10):
            name = f"{sample_index}.{id:02d}"
            dat = [
                {'method': "answer",                    'file': f"{answer_name}.{id:02d}.read..sam"},
                {'method': "vg",                        'file': f"{name}.mapped.bam"},
                # {'method': "linear",                    'file': f"{name}.index_kir_2100_raw_cons.bowtie2.bam"},
                {'method': "linear_ab",                 'file': f"{name}.index_kir_2100_ab_cons.bowtie2.bam"},
                # {'method': "full",                      'file': f"{name}.index_kir_2100_raw_full.bowtie2.bam"},
                # {'method': "ping",             'id': 0, 'file': f"data/linnil1_syn_wide.00.kir_2100_raw_full.ping.sam"},
                {'method': "hisat_merge",               'file': f"{name}.index_kir_2100_merge.mut01.bam"},
                # {'method': "hisat_merge_right",         'file': f"{name}.index_kir_2100_merge_assign1.mut01.bam"},
                # {'method': "hisat_merge_type", 'id': id, 'file': f"data/linnil1_syn_wide.{id:02d}.kir_2100_merge.mut01.hisatgenotype.sam"},
                # {'method': "hisat_raw",                 'file': f"{name}.index_kir_2100_raw.mut01.bam"},
                # {'method': "hisat_raw_type",            'file': f"{name}.index_kir_2100_raw.mut01.hisatgenotype.errcorr.sam"},
                # {'method': "hisat_raw_type_nomulti",    'file': f"{name}.index_kir_2100_raw.mut01.hisatgenotype.errcorr.no_multi.sam"},
                {'method': "hisat_2dl1s1",              'file': f"{name}.index_kir_2100_2dl1s1.mut01.bam"},
                {'method': "hisat_2dl1s1_type",         'file': f"{name}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.sam"},
                {'method': "hisat_2dl1s1_type_nomulti", 'file': f"{name}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.no_multi.sam"},
                # {'method': "hisat_splitab",             'file': f"{name}.index_kir_2100_ab.mut01.bam"},
                # {'method': "hisat_splitab_type",        'file': f"{name}.index_kir_2100_ab.mut01.hisatgenotype.errcorr.sam"},
                # {'method': "hisat_splitab_type_nomulti",'file': f"{name}.index_kir_2100_ab.mut01.hisatgenotype.errcorr.no_multi.sam"},
            ]
            for i in dat:
                i['id'] = id
            data.extend(dat)
        print(data)

        # for i in data:
        #     i.update(getStat(i['file']))
        with ProcessPoolExecutor() as executor:
            for d, result in zip(data, executor.map(getStat, map(lambda i: i['file'], data))):
                d.update(result)
        json.dump(data, open(filename + ".json", "w"))
    else:
        data = json.load(open(filename + ".json"))

    df = pd.DataFrame(data)
    ans = df[(df['method'] == 'answer')]
    ans = dict(zip(ans['id'], ans['total']))
    df['ans_total'] = df.apply(lambda i: ans[i.id], axis=1)
    df['pair_perc'] = df['pair'] / df['ans_total']
    df['secd_perc'] = df['secd'] / df['ans_total']
    # df = df[df['method'] != "hisat_merge_right"]
    # df = df[df['method'] != "full"]
    # df.loc[df['method'] == "hisat_raw", "method"] = "hisat_split"
    fig0 = px.box(df, x="method",  y="pair_perc", title="Proper Paired Ratio",
                  labels={'pair_perc': "Primary paired reads / Total reads"})
    fig1 = px.box(df, x="method",  y="secd_perc", title="Secondary Ratio",
                  labels={'secd_perc': "Secondary reads / Total reads"})
    return [fig0, fig1]


def getEachReadMappedOn(filename):
    data = defaultdict(list)
    for i in samtools("view", filename):
        if not i:
            continue
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


def customCalcMissing(total, reads, **kwargs):
    """
    Calculate the number of read cannot mapped

    Args:
      total(dict[str, int]): Nubmer(value) of reads belong to the reference(key)
      reads(dict[str, list[str]]): The gene mapped on (value) of for each read(key)
    """
    data = {gene: {'total': num, 'count': 0, 'miss': 0} for gene, num in total.items()}

    for read_name, mapping_info in reads.items():
        data[getGeneName(read_name)]['count'] += 1
        if len(list(filter(lambda i: i[1] != '*' and not (i[0] & 256), mapping_info))) < 2:
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


def customCalcAnyMapped(total, reads, **kwargs):
    """ Any of pair-reads are correctly mapped on reference """
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


def customCalcSecondary(total, reads, **kwargs):
    """ The number of secondary reads """
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


def customCalcGeneAcc(total, reads, getRenamedGeneName):
    """ The pair reads are mapped on reference """
    data = {gene: {'total': num, 'count': 0, 'primary': 0, 'secondary': 0} for gene, num in total.items()}

    for read_name, mapping_info in reads.items():
        data[getGeneName(read_name)]['count'] += 1
        # remove non-pair
        mapping_info = list(filter(lambda i: (i[0] & 2), mapping_info))

        primary = list(filter(lambda i: i[1] != '*' and not (i[0] & 256), mapping_info))
        if primary and getRenamedGeneName(primary[0][1]) == getRenamedGeneName(read_name):
            data[getGeneName(read_name)]['primary'] += 1

        second = list(filter(lambda i: i[1] != '*', mapping_info))
        if second and any(getRenamedGeneName(read[1]) == getRenamedGeneName(read_name) for read in second):
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
            'type': 'primary+secondary',
        }


def plotGenewiseMapping():
    sample_index = "data/linnil1_syn_wide.test10"
    sample_index = "data/linnil1_syn_30x"
    filename = f"{sample_index}.bam_gene_data_30x"
    answer_name = "linnil1_syn_wide/linnil1_syn_wide"
    answer_name = "linnil1_syn_30x/linnil1_syn_30x"

    if not os.path.exists(filename + ".json"):
        # Figure: Proper Paired Perecetage
        data = []
        for id in range(10):
            name = f"{sample_index}.{id:02d}"
            dat = [
                {'method': "answer",        'gene_compare_type': "",           'file': f"{answer_name}.{id:02d}.read..sam"},
                {'method': "hisat_raw",     'gene_compare_type': "ab",         'file': f"{name}.index_kir_2100_raw.mut01.bam"},
                {'method': "hisat_2dl1s1",  'gene_compare_type': "ab_2dl1s1",  'file': f"{name}.index_kir_2100_2dl1s1.mut01.bam"},
                {'method': "hisat_splitab", 'gene_compare_type': "",           'file': f"{name}.index_kir_2100_ab.mut01.bam"},
                {'method': "linear",        'gene_compare_type': "ab",         'file': f"{name}.index_kir_2100_raw_cons.bowtie2.bam"},
                {'method': "linear_ab",     'gene_compare_type': "",           'file': f"{name}.index_kir_2100_ab_cons.bowtie2.bam"},
                # {'method': "full",   'id': 0, 'file': "data/linnil1_syn_wide.00.kir_2100_raw_full.bowtie2.sam"},
            ]
            for i in dat:
                i['id'] = id
            data.extend(dat)

        # for i in data:
        #     i.update(getEachReadMappedOn(i['file']))
        with ProcessPoolExecutor() as executor:
            for d, result in zip(data, executor.map(getEachReadMappedOn, map(lambda i: i['file'], data))):
                d.update(result)

        # print(data)
        json.dump(data, open(filename + ".json", "w"))
    else:
        print('Reading')
        data = json.load(open(filename + ".json"))

    def reColor(arr1, arr2):
        """
        Assign color for arr1_i x arr2_j

        Args:
          arr1(list[str]): ...
          arr2(list[str]): ...
        Returns:
          group_name(list[str]):  list of arr1_i x arr2_j
          name_to_color(dict[str, str]): key=arr1_i-arr2_j value=color
        """
        group_name = arr1 + " - " + arr2
        colors = (px.colors.qualitative.Dark2, px.colors.qualitative.Set2)

        name_to_color = {}
        for i, m in enumerate(set(arr1)):
            for j, t in enumerate(set(arr2)):
                name_to_color[m + " - " + t] = colors[j][i]
        return group_name, name_to_color

    def findAnswerSample(id):
        """ Find the sample dict with same id """
        for sample in data:
            if sample['method'] == "answer" and sample['id'] == id:
                return sample
        raise ValueError

    def runAllFiles(data, stat_func):
        """
        Calculate the stat by stat_func for each sample(data)

        Args:
          data(dict[str, Any]): The data from
        """
        custom_stats = []
        for sample in data:
            if sample['method'] == "answer":
                continue
            print(sample['method'], sample['id'])
            reads_ans = findAnswerSample(sample['id'])['mapping']
            reads_file = sample['mapping']
            gene_ans_total = Counter(map(lambda name: getGeneName(name), reads_ans.keys()))

            if sample['gene_compare_type'] == "":
                stats = stat_func(gene_ans_total, reads_file, getRenamedGeneName=getGeneName)
            elif sample['gene_compare_type'] == "ab":
                stats = stat_func(gene_ans_total, reads_file, getRenamedGeneName=lambda i: i[:7])
            elif sample['gene_compare_type'] == "ab_2dl1s1":
                def compare2dl1s1(name):
                    name = getGeneName(name)
                    return {
                        "KIR2DL5A": "KIR2DL5",
                        "KIR2DL5B": "KIR2DL5",
                        "KIR2DL1": "KIR2DL1S1",
                        "KIR2DS1": "KIR2DL1S1",
                    }.get(name, name)
                stats = stat_func(gene_ans_total, reads_file, getRenamedGeneName=compare2dl1s1)

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
    df = runAllFiles(data, customCalcMissing)
    gene_order = {'gene': sorted(set(df['gene']))}
    figs.append(
        px.box(df, x="gene", y="miss_perc", color="method", category_orders=gene_order)
          .update_layout(title="Missed reads belong to (Proportion)",
                         yaxis_title="Missed_reads / total_read_counts in that gene",
                         yaxis_tickformat='e')
    )
    figs.append(
        px.box(df, x="method", y="miss_perc", color="method", category_orders=gene_order)
          .update_layout(title="Missed reads Per Method (Proportion)",
                         yaxis_title="Missed_reads / total_read_counts in that gene",
                         yaxis_tickformat='e')
    )

    # Any proper reads
    df = runAllFiles(data, customCalcAnyMapped)
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
    df = runAllFiles(data, customCalcSecondary)
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
    df = runAllFiles(data, customCalcGeneAcc)
    gene_order = {'gene': sorted(set(df['gene']))}
    group, color = reColor(df["method"], df["type"])
    figs.append(
        px.box(df, x="gene", y="acc", color=group,
               category_orders=gene_order,
               color_discrete_map=color,
               ).update_layout(
                   title="Gene-level accuracy(Primary)",
                   yaxis_title="Correct read / total read (%)")
    )
    figs.append(
        px.box(df, x="method", y="acc", color="type", category_orders=gene_order)
          .update_layout(title="Gene-level accuracy Per Method (Proportion)",
                         yaxis_title="Correct read / total read (%)")
    )
    return figs


def extractTargetReadFromBam(bamfile, mapped_ref, wanted_gene):
    suffix = f".extract.{mapped_ref.split('*')[0]}.{wanted_gene}"
    outputfile = os.path.splitext(bamfile)[0] + suffix

    # read
    data = []
    for i in samtools("view", f"-h {bamfile}"):
        if not i:
            continue
        if i.startswith("@"):
            data.append(i)
        elif i.split("\t")[2] == mapped_ref and wanted_gene in i.split('*')[0]:
            data.append(i)
            continue
    with open(outputfile + ".sam", "w") as f:
        f.writelines(map(lambda i: i + "\n", data))
    samtobam(outputfile)
    print("Save to", outputfile + ".bam")
    return suffix


if __name__ == '__main__':
    """
    extractTargetReadFromBam("data/linnil1_syn_wide.test10.00.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.no_multi.sam",
                             mapped_ref="KIR2DL1S1*BACKBONE",
                             wanted_gene="KIR2DL2")
    extractTargetReadFromBam("data/linnil1_syn_wide.test10.00.index_kir_2100_raw.mut01.hisatgenotype.errcorr.no_multi.sam",
                             mapped_ref="KIR2DL1*BACKBONE",
                             wanted_gene="KIR2DS1")
    """
    figs = []
    # figs.extend(plotBamStat())
    # figs.extend(plotGenewiseMapping())
    app = Dash(__name__)
    app.layout = html.Div([dcc.Graph(figure=f) for f in figs])
    app.run_server(debug=True, port=8052)
