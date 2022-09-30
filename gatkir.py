# the pipeline is modified from https://github.com/saorisakaue/KIR_project
import os
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from collections import Counter

from namepipe import nt, NameTask
from kg_utils import (
    runShell,
    runDocker,
    threads,
    samtobam,
)
from kg_eval import readAnswerAllele, compareCohort


def getGeneName(s):
    return s.split("*")[0]


@nt
def linkSamples(input_name, data_folder):
    # 1 -> 1
    name = input_name.split('/')[0]
    output_name = os.path.join(data_folder, name + ".{}")
    new_name = output_name.format(input_name.template_args[0])
    if Path(new_name + ".read.1.fq").exists():
        return output_name
    runShell(f"ln -s ../{input_name}.1.fq {new_name}.read.1.fq")
    runShell(f"ln -s ../{input_name}.2.fq {new_name}.read.2.fq")
    return output_name


@nt
def setupGATKIR(input_name):
    folder = "gatkir"
    if Path(folder).exists():
        return folder
    runShell(f"git clone https://github.com/saorisakaue/KIR_project {folder}")
    runShell('ln -s "KIR2DL5A;KIR2DL5B.difpos.all.txt" gatkir/data/KIR2DL5AB.difpos.all.txt')
    runShell('ln -s "KIR2DS3;KIR2DS5.difpos.all.txt" gatkir/data/KIR2DS35.difpos.all.txt')
    return folder


@nt
def bwa(input_name, index):
    id = input_name.template_args[0]
    output_name = input_name + ".bwa"
    if Path(output_name + ".bam").exists():
        return output_name

    rg = "@RG\\tID:" + id + "\\tSM: " + id
    runDocker("bwa", f""" \
        bwa mem -t {threads} \
            {index}/REF/KIR_seq_ref \
            -R "{rg}" \
            {input_name}.read.1.fq \
            {input_name}.read.2.fq \
            -o {output_name}.sam
    """)
    samtobam(output_name)
    return output_name


@nt
def addGroup(input_name):
    output_name = input_name + ".rg"
    id = input_name.template_args[0]
    if Path(output_name + ".bam").exists():
        return output_name
    runDocker("picard", f""" \
        picard AddOrReplaceReadGroups \
        I={input_name}.bam \
        O={output_name}.bam \
        RGLB={id} RGPL=ILLUMINA RGPU={id} RGSM={id} RGID={id} \
        VALIDATION_STRINGENCY=LENIENT
    """)
    return output_name


@nt
def markDupliate(input_name):
    output_name = input_name + ".md"
    if Path(output_name + ".bam").exists():
        return output_name
    runDocker("picard", f""" \
        picard MarkDuplicates \
        I={input_name}.bam \
        O={output_name}.bam \
        ASSUME_SORTED=false \
        REMOVE_DUPLICATES=false \
        CREATE_INDEX=True \
        VALIDATION_STRINGENCY=LENIENT \
        M={output_name}.metrics
    """)
    return output_name


@nt
def analysisTK(input_name, index):
    output_name = input_name + ".coverage"
    if Path(output_name + ".vcf").exists():
        return output_name
    runDocker("gatk3", f""" \
        java -Xmx40g -Xms40g -jar /usr/GenomeAnalysisTK.jar \
        -T DiagnoseTargets \
        -I {input_name}.bam \
        -o {output_name}.vcf \
        -R {index}/REF/KIR_seq_ref.fasta \
        -L {index}/REF/KIR_seq_ref.intervals
    """)
    return output_name


@nt
def getCoverage(input_name):
    output_name = input_name + ".depth_per_gene"
    if Path(output_name + ".csv").exists():
        return output_name
    df = []
    for line in open(input_name + ".vcf"):
        if not line or line.startswith("#"):
            continue
        # KIR2DL1	1	.	C	<DT>	.	PASS	END=159;IDP=17.63;IGC=0.604	IDP:LL:ZL	17.63:0:0
        line = line.split("\t")
        info = dict(map(lambda i: i.split("="), line[7].split(';')))
        gene = line[0]
        length = float(info["END"]) - float(line[1])
        df.append({
            'gene': gene,
            'depth': float(info["IDP"]),
            'length': length
        })

    # linnil1: average depth with weighted by length
    df = pd.DataFrame(df)
    df = df.groupby("gene").apply(lambda i: np.average(i.depth, weights=i.length))
    df = df.reset_index()
    print(df)
    df.to_csv(output_name + ".csv", index=False, header=None)
    return output_name


def getThreshold(cov, name):
    # 02_determine_copy_number.py
    colors = ['r','g','b','c','m','y']

    thres_dic = {}
    for i in range(len(cov)):
        gene = cov.index[i]
        a = np.array(cov.iloc[i,:]).reshape(-1,1)
        kde = KernelDensity(kernel='gaussian', bandwidth=0.075).fit(a)
        s = np.linspace(min(a)-0.05,max(a)+0.05)
        e = kde.score_samples(s.reshape(-1,1))
        mi, ma = argrelextrema(e,np.less)[0], argrelextrema(e,np.greater)[0]
        n_thress = len(mi)
        n_bin = len(ma)
        thres_dic[gene] = s[mi]
        if n_thress > 0:
            plt.figure()
            plt.plot(s[:mi[0]+1],e[:mi[0]+1], colors[0])
            for j in range(1,n_bin-1):
                plt.plot(s[mi[j-1]:mi[j]+1],e[mi[j-1]:mi[j]+1], colors[j])
            plt.plot(s[mi[n_bin-2]:],e[mi[n_bin-2]:],colors[(n_bin-1) % len(colors)])
            plt.savefig(f"{name}.{gene}.png")
            print(f"Save {name}.{gene}.png")
        else:
            print(gene + ' had zero threshold')
    return thres_dic


def getPloidy(cov, thres_dic):
    # 02_determine_copy_number.py
    genelist = ['KIR3DS1', 'KIR3DL1', 'KIR2DS4', 'KIR2DS35', 'KIR2DS2',
                'KIR2DS1', 'KIR2DP1', 'KIR2DL5AB', 'KIR2DL3', 'KIR2DL2',
                'KIR2DL1', 'KIR3DL3', 'KIR3DL2', 'KIR2DL4']
    tmp = np.zeros((len(genelist),len(cov.columns)),dtype = "int")
    copy = pd.DataFrame(tmp)
    copy.columns = cov.columns
    copy.index = genelist
    for gene in genelist:
        cut_thres = np.hstack(([0], np.ravel(thres_dic[gene]), [4]))
        ratio = cov.loc[gene, :]
        copy.loc[gene, :] = np.array(pd.cut(ratio, cut_thres, labels=False))
    return copy


@nt
def ploidyEstimate(input_name):
    names = input_name.get_input_names()
    output_name = input_name.replace_wildcard("_merge_depth") + ".ploidy"

    if Path(output_name + ".csv").exists():
        return output_name

    # merge df
    dfs = []
    for name in names:
        df = pd.read_csv(name + ".csv", header=None, index_col=0)
        df = df.rename(index={"KIR2DL5A;KIR2DL5B": "KIR2DL5AB", "KIR2DS3;KIR2DS5": "KIR2DS35"})
        df.columns = [name.template_args[0]]
        dfs.append(df)
    df = pd.concat(dfs, axis=1)

    # normalize by 3DL3
    df = df / df.loc["KIR3DL3", :]
    print(df)

    # 02_determine_copy_number
    thres_dic = getThreshold(df, input_name.replace_wildcard("_merge_depth"))
    ploidy = getPloidy(df, thres_dic)
    ploidy.loc["KIR3DL3", :] = 2
    ploidy = ploidy.fillna(0).astype(int)
    print(thres_dic, ploidy)

    # save
    ploidy.to_csv(output_name + ".csv")
    return output_name


@nt
def ploidyExploded(input_name, ploidy_name):
    output_name = input_name + ".ploidy_" + ploidy_name.split('/')[-1].replace(".", "_") + ".{}"
    if Path(output_name.format("KIR3DL3") + ".txt").exists():
        return output_name

    ploidy = pd.read_csv(ploidy_name + ".csv", index_col=0)
    sample_ploidy = ploidy[input_name.template_args[0]]
    print(sample_ploidy)

    for gene, ploidy in sample_ploidy.iteritems():
        if not ploidy:
            continue
        if gene == "KIR3DP1":
            continue
        runShell(f"ln -s ../{input_name}.bam {output_name.format(gene)}.bam")
        runShell(f"ln -s ../{input_name}.bai {output_name.format(gene)}.bai")
        open(f"{output_name.format(gene)}.txt", "w").write(str(ploidy))
    return output_name


@nt
def haplotypeCaller(input_name, index):
    output_name = input_name + ".hc"
    ploidy = int(open(input_name + ".txt").read())
    gene = input_name.template_args[1]
    # Note: this is gvcf
    if Path(output_name + ".g.vcf.gz").exists():
        return output_name

    print(input_name, ploidy)
    runDocker("gatk3", f""" \
        java -Xmx40g -Xms40g -jar /usr/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -I {input_name}.bam \
        -o {output_name}.g.vcf.gz \
        -nct 2 \
        -ploidy {ploidy} \
        -R {index}/REF/KIR_seq_ref.fasta \
        -L {index}/REF/{gene}.intervals \
        --emitRefConfidence GVCF
    """)
    return output_name


@nt
def jointGenotype(input_name, index):
    output_name = input_name.replace_wildcard("_jg")

    files = [name + ".g.vcf.gz" for name in input_name.get_input_names()]

    if Path(output_name + ".g.vcf.gz").exists():
        return output_name

    runDocker("gatk3", f""" \
        java -Xmx40g -Xms40g -jar /usr/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -R {index}/REF/KIR_seq_ref.fasta \
        -allSites \
        -o {output_name}.g.vcf.gz \
        {" ".join(map(lambda i: "--variant " + i, files))}
    """)
    return output_name


@nt
def separteVCFSample(input_name, genotype):
    output_name = input_name + ".from_" + genotype.split('/')[-1].replace('.', "_")
    if Path(output_name + ".g.vcf.gz").exists():
        return output_name

    runDocker("bcftools", f""" \
        bcftools view {genotype}.g.vcf.gz -s {input_name.template_args[0]} -o {output_name}.g.vcf.gz
    """)
    return output_name


@nt
def deepVariant(input_name, index):
    output_name = input_name + ".dv"
    if Path(output_name + ".g.vcf.gz").exists():
        return output_name

    runDocker("deepvariant", f""" \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref {index}/REF/KIR_seq_ref.fasta \
        --reads {input_name}.bam \
        --output_vcf={output_name}.vcf.gz \
        --output_gvcf={output_name}.g.vcf.gz \
     """)
    return output_name


@nt
def vcfNorm(input_name, index):
    output_name = input_name + ".norm"
    if Path(output_name + ".g.vcf.gz").exists():
        return output_name

    runDocker("bcftools", f""" \
        bcftools norm \
        --fasta-ref {index}/REF/KIR_seq_ref.fasta \
        {input_name}.g.vcf.gz -o {output_name}.g.vcf.gz
    """)
    return output_name


@nt
def calling(input_name, index):
    output_name = input_name + ".call"
    if Path(output_name + ".KIR3DL3.alleles.tsv").exists():
        return output_name + ".{}"

    id = input_name.template_args[0]
    GENE_REF = {
        'KIR2DL1': 'KIR2DL1_001',
        'KIR2DL2': 'KIR2DL2_0010101',
        'KIR2DL3': 'KIR2DL3_0010101',
        'KIR2DL5AB': 'KIR2DL5A_0010101',
        'KIR2DS1': 'KIR2DS1_001',
        'KIR2DS2': 'KIR2DS2_0010101',
        'KIR2DS35': 'KIR2DS3_00101',
        'KIR2DS4': 'KIR2DS4_0010101',
        'KIR3DL1': 'KIR3DL1_0010101',
        'KIR3DL2': 'KIR3DL2_0010101',
        'KIR3DL3': 'KIR3DL3_00101',
        'KIR3DS1': 'KIR3DS1_010',
        'KIR2DL4': 'KIR2DL4_00101'
    }

    for gene in GENE_REF:
        runShell(f"""
            python gatkir_call.py \
            {input_name}.g.vcf.gz {index}/data/{gene}.difpos.all.txt \
            {output_name}.{gene}.dosage.tsv {output_name}.{gene}.reference.tsv {output_name}.{gene}.alleles.tsv \
            {gene} {id}
        """)
    return output_name + ".{}"


@nt
def mergeCall(input_name):
    files = [i + ".alleles.tsv" for i in input_name.get_input_names()]
    output_name = input_name.replace_wildcard("_merge")
    if Path(output_name + ".tsv").exists():
        return output_name
    runShell(f"cat {' '.join(files)} > {output_name}.tsv")
    return output_name


def readGatkirAlleles(filename, select_all=False):
    raw_df = pd.read_csv(filename, header=None, sep="\t")
    raw_df.columns = ["id", "gene", "alleles", "type"]
    alleles = []
    for i in raw_df.itertuples():
        if i.type == "known":
            alleles_text = i.alleles.replace("_", "*")
            possible_set = alleles_text.split("-or-")
        elif i.type == "potentially_novel":
            alleles_text = i.alleles.replace("Close_to_", "").replace("_", "*").split("[")[0]
            possible_set = alleles_text.split("-OR-")
        else:
            raise "type not found"

        if select_all:
            alleles.extend([k for i in possible_set for j in i.split("/") for k in j.split('-')])
        else:
            alleles.extend([i.split('-')[0] for i in possible_set[0].split("/")])
    return alleles


@nt
def mergeAlleles(input_name, answer, select_all=False):
    if select_all:
        output_name = input_name.replace_wildcard("_merge_called_full")
    else:
        output_name = input_name.replace_wildcard("_merge_called")
    called_alleles_dict = {}

    for name in input_name.get_input_names():
        alleles = readGatkirAlleles(name + ".tsv", select_all=select_all)
        id = name.template_args[0]
        called_alleles_dict[id] = alleles
        print(id, alleles)

    predict_list = []
    for id, alleles in called_alleles_dict.items():
        predict_list.append({
            'id': id,
            'alleles': '_'.join(alleles),
            'name': input_name.format(id),
        })
    df = pd.DataFrame(predict_list)
    df.to_csv(f"{output_name}.tsv", index=False, sep="\t")
    print(df)

    ans = readAnswerAllele(f"{answer}/{answer}.summary.csv")
    compareCohort(ans, called_alleles_dict)
    return output_name


@nt
def answerPloidy(folder_name, answer):
    output_name = f"{folder_name}/{answer}_answer_cn"
    if Path(output_name + ".csv").exists():
        return output_name

    ans = readAnswerAllele(f"{answer}/{answer}.summary.csv")
    ids = []
    dfs = []
    def renameGene(i):
        return {
            'KIR2DS3': "KIR2DS35",
            'KIR2DS5': "KIR2DS35",
            'KIR2DL5A': "KIR2DL5AB",
            'KIR2DL5B': "KIR2DL5AB",
        }.get(i, i)

    for id, alleles in ans.items():
        ids.append(id)
        dfs.append(Counter(map(renameGene, map(getGeneName, alleles))))

    df = pd.DataFrame(dfs)
    df.index = ids
    df = df.T
    df = df.fillna(0).astype(int)
    df.to_csv(output_name + ".csv")
    print(df)

    return output_name


def testThreshold():
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from sklearn.neighbors import KernelDensity
    from scipy.signal import argrelextrema
    from graphkir.gk_cn_model import KDEcut
    from graphkir.gk_plot import showPlot

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
    kde = KDEcut()
    kde.fit(ori_data)
    x = np.linspace(0, 1.1, kde.points)
    y = kde.kde.score_samples(x[:, None])
    fig.add_trace(go.Scatter(x=x * kde.x_max / norm_value, y=y / np.abs(y.min()),  name=f"our_method", line_color="blue"))
    for cn, m in enumerate(kde.local_min):
        fig.add_vline(x=m * kde.x_max / norm_value, line_width=2, line_dash="dash", line_color="blue", annotation_text=f"CN={cn}", annotation_yshift=-20, annotation_font_color="blue")

    # Bad paramerts: bandwidth
    data = np.array(ori_data)[:, None] / norm_value
    x = np.linspace(data.min() - 0.05, data.max() + 0.05)
    kde  = KernelDensity(kernel='gaussian', bandwidth=0.025).fit(data)
    y  = kde.score_samples(x[:, None])
    mi, ma   = argrelextrema(y, np.less)[0], argrelextrema(y, np.greater)[0]
    fig.add_trace(go.Scatter(x=x, y=y / np.abs(y.min()),  name=f"Bad parameters", line_color="red"))
    for cn, m in enumerate(mi):
        fig.add_vline(x=x[m], line_width=2, line_dash="dash", line_color="red", annotation_text=f"CN={cn}", annotation_font_color="red")

    showPlot([fig])


def testShowIdea():
    from ping import readPingLocusCount
    import plotly.express as px
    from graphkir.gk_cn_model import KDEcut
    from graphkir.gk_plot import showPlot

    folder = "data3/ping_linnil1_syn_30x_seed87.result"
    df_ping = readPingLocusCount(f"{folder}/locusRatioFrame.csv")
    threshold_df = pd.read_csv(f"{folder}/manualCopyThresholds.csv")
    threshold_df.columns = ["gene", *threshold_df.columns[1:]]
    df = pd.concat([df_ping])
    print(df)
    print(threshold_df)
    df = df.melt(["id", "method"], var_name="gene")
    df = df.sort_values(["method", "value"], ascending=[False, True])  # PING value is ascending

    # plot
    figs = []
    for gene in sorted(set(df["gene"])):
        df_gene = df[df["gene"] == gene]
        fig = px.scatter(df_gene, x="id", y="value", title=gene)
        fig.update_layout(yaxis_title=f"{gene}/KIR3DL3 ratio",
                          xaxis_title="Sample ID")
        figs.append(fig)
        kde = KDEcut()
        kde.bandwidth = 0.08
        kde.fit(df_gene['value'])
        figs.extend(kde.plot())
        figs[-1].update_layout(title=f"{gene}'s copy number estimation by KDE")

    showPlot(figs)


if __name__ == "__main__":
    # testThreshold()
    # testShowIdea()
    # exit()
    data_folder = "data3"
    answer = "linnil1_syn_30x_seed87"
    Path(data_folder).mkdir(exist_ok=True)
    index = None >> setupGATKIR
    samples = answer + "/" + answer + ".{}.read" >> linkSamples.set_args(data_folder)
    mapping = samples >> bwa.set_args(index=index) \
                      >> addGroup \
                      >> markDupliate
    # switch this
    ploidy = mapping >> analysisTK.set_args(index=index) \
                     >> getCoverage \
                     >> ploidyEstimate.set_depended(-1)
    ploidy = data_folder >> answerPloidy.set_args(answer=answer)

    genotype_joint = mapping >> ploidyExploded.set_args(ploidy_name=str(ploidy)) >> \
                                haplotypeCaller.set_args(index=index) >> \
                                jointGenotype.set_args(index=index).set_depended([0, 1])
    gatk = mapping >> separteVCFSample.set_args(genotype=str(genotype_joint)) \
                   >> vcfNorm.set_args(index=index)
    gatk >> calling.set_args(index=index) \
         >> mergeCall.set_depended(-1) \
         >> mergeAlleles.set_args(answer=answer, select_all=True).set_depended(-1)

    # TODO
    # deep_variant = mapping >> deepVariant.set_args(index=index) \
    #                        >> vcfNorm.set_args(index=index)
    # How to use deep_variant in gatkir_call
