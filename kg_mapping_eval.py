from itertools import chain
import os
import json
import subprocess
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from dash import Dash, dcc, html, Input, Output
import plotly.graph_objects as go
import plotly.express as px

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
    answer = "linnil1_syn_30x_seed87"
    filename = f"data/{answer}.bam_stat"
    # filename = f"{sample_index}.bam_stat_vg"

    if not os.path.exists(filename + ".json"):
        # Figure: Proper Paired Perecetage
        id = 0
        data = []
        for id in range(10):
            name = f"data/{answer}.{id:02d}"
            dat = [
                {'method': "answer",                       'file': f"{answer}/{answer}.{id:02d}.read..sam"},
                {'method': "linear_ab",                    'file': f"{name}.index_kir_2100_raw_cons.bowtie2.bam"},
                {'method': "linear",                       'file': f"{name}.index_kir_2100_ab_cons.bowtie2.bam"},
                {'method': "linear_ab_2dl1s1",             'file': f"{name}.index_kir_2100_2dl1s1_cons.bowtie2.bam"},
                {'method': "linear_full",                  'file': f"{name}.index_kir_2100_raw_full.bowtie2.bam"},
                {'method': "bwa_ab",                       'file': f"{name}.index_kir_2100_raw_cons_bwa.bwa.bam"},
                {'method': "bwa_ab_2dl1s1",                'file': f"{name}.index_kir_2100_2dl1s1_cons_bwa.bwa.bam"},
                {'method': "bwa",                          'file': f"{name}.index_kir_2100_ab_cons_bwa.bwa.bam"},
                {'method': "bwa_merge",                    'file': f"{name}.index_kir_2100_merge_cons_bwa.bwa.bam"},
                {'method': "hisat_merge",                  'file': f"{name}.index_kir_2100_merge.mut01.bam"},
                {'method': "hisat_merge_type",             'file': f"{name}.index_kir_2100_merge.mut01.hisatgenotype.errcorr.bam"},
                {'method': "hisat_merge_type_nomulti",     'file': f"{name}.index_kir_2100_merge.mut01.hisatgenotype.errcorr.no_multi.bam"},
                {'method': "hisat_ab",                     'file': f"{name}.index_kir_2100_raw.mut01.bam"},
                {'method': "hisat_ab_type",                'file': f"{name}.index_kir_2100_raw.mut01.hisatgenotype.errcorr.bam"},
                {'method': "hisat_ab_type_nomulti",        'file': f"{name}.index_kir_2100_raw.mut01.hisatgenotype.errcorr.no_multi.bam"},
                {'method': "hisat_ab_2dl1s1",              'file': f"{name}.index_kir_2100_2dl1s1.mut01.bam"},
                {'method': "hisat_ab_2dl1s1_type",         'file': f"{name}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.bam"},
                {'method': "hisat_ab_2dl1s1_type_nomulti", 'file': f"{name}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.no_multi.bam"},
                {'method': "hisat",                        'file': f"{name}.index_kir_2100_ab.mut01.bam"},
                {'method': "hisat_type",                   'file': f"{name}.index_kir_2100_ab.mut01.hisatgenotype.errcorr.bam"},
                {'method': "hisat_type_nomulti",           'file': f"{name}.index_kir_2100_ab.mut01.hisatgenotype.errcorr.no_multi.bam"},
                {'method': "hisat_271-ab-2dl1s1",          'file': f"{name}.index_kir_271_2dl1s1.mut01.bam"},
                {'method': "hisat_290-ab-2dl1s1",          'file': f"{name}.index_kir_290_2dl1s1.mut01.bam"},
                {'method': "vg_ab",                        'file': f"data3/{answer}.{id:02d}.data3_kir_2100_raw.vgindex.bam"},
                {'method': "vg_merge",                     'file': f"data3/{answer}.{id:02d}.data3_kir_2100_merge.vgindex.bam"},
                {'method': "ping_allele_setup",            'file': f"data3/ping_{answer}.result/allele_setup_files/{answer}.{id:02d}.read..bam"},
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
    # sample_index = "data/linnil1_syn_wide.test10"
    answer = "linnil1_syn_30x_seed87"
    sample_index = f"data/{answer}"
    # filename = f"{sample_index}.bam_gene_data_30x"
    filename = "123"

    if not os.path.exists(filename + ".json"):
        # Figure: Proper Paired Perecetage
        data = []
        for id in range(10):
            name = f"{sample_index}.{id:02d}"
            dat = [
                {'method': "answer",        'gene_compare_type': "",
                 'file': f"{answer}/{answer}.{id:02d}.read..sam"},
                {'method': "hisat_ab",                     'gene_compare_type': "ab",        'file': f"{name}.index_kir_2100_raw.mut01.bam"},
                {'method': "hisat_ab_type",                'gene_compare_type': "ab",        'file': f"{name}.index_kir_2100_raw.mut01.hisatgenotype.errcorr.bam"},
                {'method': "hisat_ab_type_nomulti",        'gene_compare_type': "ab",        'file': f"{name}.index_kir_2100_raw.mut01.hisatgenotype.errcorr.no_multi.bam"},
                {'method': "hisat",                        'gene_compare_type': "",          'file': f"{name}.index_kir_2100_ab.mut01.bam"},
                {'method': "hisat_type",                   'gene_compare_type': "",          'file': f"{name}.index_kir_2100_ab.mut01.hisatgenotype.errcorr.bam"},
                {'method': "hisat_type_nomulti",           'gene_compare_type': "",          'file': f"{name}.index_kir_2100_ab.mut01.hisatgenotype.errcorr.no_multi.bam"},
                {'method': "hisat_ab_2dl1s1",              'gene_compare_type': "ab_2dl1s1", 'file': f"{name}.index_kir_2100_2dl1s1.mut01.bam"},
                {'method': "hisat_ab_2dl1s1_type",         'gene_compare_type': "ab_2dl1s1", 'file': f"{name}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.bam"},
                {'method': "hisat_ab_2dl1s1_type_nomulti", 'gene_compare_type': "ab_2dl1s1", 'file': f"{name}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.no_multi.bam"},
                {'method': "linear_ab",                    'gene_compare_type': "ab",         'file': f"{name}.index_kir_2100_raw_cons.bowtie2.bam"},
                {'method': "linear",                       'gene_compare_type': "",           'file': f"{name}.index_kir_2100_ab_cons.bowtie2.bam"},
                {'method': "linear_ab_2dl1s1",             'gene_compare_type': "ab_2dl1s1",  'file': f"{name}.index_kir_2100_2dl1s1_cons.bowtie2.bam"},
                {'method': "bwa_ab",                       'gene_compare_type': "ab",         'file': f"{name}.index_kir_2100_raw_cons_bwa.bwa.bam"},
                {'method': "bwa",                          'gene_compare_type': "",           'file': f"{name}.index_kir_2100_ab_cons_bwa.bwa.bam"},
                {'method': "bwa_ab_2dl1s1",                'gene_compare_type': "ab_2dl1s1",  'file': f"{name}.index_kir_2100_2dl1s1_cons_bwa.bwa.bam"},
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
                name_to_color[m + " - " + t] = colors[j % len(colors[0])][i % len(colors[1])]
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


def plotGeneFromToBar(df, x, y):
    order = {x: sorted(set(df[x])), y: sorted(set(df[y]))}
    df['txt'] = df[y] + ["*" if i else "" for i in df["is_multi"]]
    return px.bar(df, x=x, y="value", color=y, text='txt',
                  category_orders=order, color_discrete_sequence=px.colors.qualitative.Light24)


def plotGeneFromTo(df, title, x="from"):
    df = df.rename(columns={"size": "value"})
    figs = [
        plotGeneFromToBar(df, x, "to"),
        plotGeneFromToBar(df[df['is_multi'] == False], x, "to"),
        # reverse
        plotGeneFromToBar(df, "to", x),
        plotGeneFromToBar(df[df['is_multi'] == False], "to", x),
    ]
    for fig in figs:
        fig.update_layout(title=title)

    return figs


def isGeneCorrect(a, b):
    if a == "KIR2DL5":
        # b = "KIR2DL5A" or "KIR2DL5B"
        return a in b
    elif a == "KIR2DL1S1":
        return b in ["KIR2DL1", "KIR2DS1"]
    return a == b


def plotGeneMappedAcc():
    from kg_typping import getNH

    # sample_index = "data/linnil1_syn_wide.test10"
    answer = "linnil1_syn_30x_seed87"
    sample_index = f"data/{answer}"
    # filename = f"{sample_index}.bam_gene_data_30x"

    # Figure: Proper Paired Perecetage
    data = []
    for id in range(10):
        name = f"{sample_index}.{id:02d}"
        dat = [
            {'method': "hisat",           'file': f"{name}.index_kir_2100_raw.mut01.hisatgenotype.errcorr.id_only.json"},
            {'method': "hisat_ab_2dl1s1", 'file': f"{name}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.id_only.json"},
            {'method': "hisat_ab",        'file': f"{name}.index_kir_2100_ab.mut01.hisatgenotype.errcorr.id_only.json"},
        ]
        for i in dat:
            i['id'] = id
        data.extend(dat)

    figs = []
    df_multi_acc = []
    df_nomulti_acc = []
    for i in data:
        print(i)
        data = json.load(open(i['file']))
        df_ori = []
        for backbone, reads in data['reads'].items():
            gene = getGeneName(backbone)
            for read in reads:
                ans_gene = getGeneName(read['l_sam'])
                df_ori.append({
                    'from_allele': getAlleleName(read['l_sam']),
                    'from': ans_gene,
                    'to': gene,
                    'correct': isGeneCorrect(gene, ans_gene),
                    'is_multi': getNH(read['l_sam']) > 1,
                    'multi': getNH(read['l_sam']),
                })

        df_ori = pd.DataFrame(df_ori)

        # gene-level
        """
        df = df_ori.groupby(["from", "to", "is_multi", "correct"], as_index=False).size()
        figs.extend(plotGeneFromTo(df, i['file'], "from"))

        # allele-level
        df = df_ori.groupby(["from_allele", "to", "is_multi", "correct"], as_index=False).size()
        figs.extend(plotGeneFromTo(df, i['file'], "from_allele"))
        """

        # Ratio of correct gene
        df = df_ori.groupby(["to", "correct"]).size()
        df.name = "acc"
        df = df / df.groupby("to").sum(0)
        df = df.reset_index()
        df = df[df['correct'] == True]
        df['method'] = i['method']
        df_multi_acc.append(df)

        # Ratio of correct gene without mutliple alignment
        df = df_ori[df_ori['is_multi'] == False].groupby(["to", "correct"]).size()
        df.name = "acc"
        df = df / df.groupby("to").sum(0)
        df = df.reset_index()
        df = df[df['correct'] == True]
        df['method'] = i['method']
        df_nomulti_acc.append(df)

    df1 = pd.concat(df_multi_acc)
    df1['multi'] = True
    df2 = pd.concat(df_nomulti_acc)
    df2['multi'] = False
    df = pd.concat([df1, df2])
    df['method-multi'] = df['method'] + "-" + df['multi'].apply(lambda i: "multi" if i else "nomulti")
    order = {
        'method-multi': sorted(set(df['method-multi']), key=lambda i: (i.split('-')[0], i.split('-')[1] == "nomulti")),
        'to': sorted(set(df['to'])),
    }
    color = list(chain.from_iterable(zip(px.colors.qualitative.Pastel2, px.colors.qualitative.Dark2)))
    figs.append(px.box(df, x="to", y="acc", color="method-multi", category_orders=order, color_discrete_sequence=color)
                  .update_layout(yaxis_title="Gene-level correct reads / all reads on the gene"))

    return figs


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
    figs.extend(plotBamStat())
    # figs.extend(plotGenewiseMapping())
    # figs.extend(plotGeneMappedAcc())
    app = Dash(__name__)
    app.layout = html.Div([dcc.Graph(figure=f) for f in figs])
    app.run_server(debug=True, port=8052)
