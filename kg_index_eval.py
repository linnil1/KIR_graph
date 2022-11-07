import os
import json
from pprint import pprint
from typing import Iterable
from itertools import chain
from functools import reduce
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dash import Dash, dcc, html, Input, Output
from pyhlamsa import msaio, Genemsa

from graphkir.kir_msa import readFromMSAs, muscle
from graphkir.utils import runShell, samtobam, getGeneName, threads, getAlleleField
from kg_utils import runDocker


def calculate400bpReadMatchness(a1: str, a2: str) -> pd.DataFrame:
    """ deprecated: Mark if the 400-bp segment in a1 can be mapped on a2 """
    # read
    index = "index/kir_2100_merge.save"
    output_name = f"{index}.400bp_match.{a1}_{a2}.csv"
    if os.path.exists(output_name):
        return pd.read_csv(output_name)
    msa = msaio.load_msa(f"{index}.KIR.fa", f"{index}.KIR.json")

    # separate
    msa = msa.select_allele(f"{a1}.*|{a2}.*").shrink().reset_index().sort_name()
    msa_a1 = [(i[0], i[1].replace("-", "")) for i in msa.select_allele(f"{a1}.*").alleles.items()]
    msa_a2 = [(i[0], i[1].replace("-", "")) for i in msa.select_allele(f"{a2}.*").alleles.items()]

    def isIn(items: Iterable[tuple[str, str]], seq: str) -> bool | str:
        for i in items:
            if seq in i[1]:
                return i[0]
        return False

    # main
    df = pd.DataFrame()
    for name, seq in msa_a1:
        df.loc[name, "name"] = name
        for pos in range(0, len(seq) - 100, 100):
            allele_finded = isIn(msa_a2, seq[pos: pos+400])
            if allele_finded:
                print(name, pos, allele_finded)
            df.loc[name, pos] = not (not allele_finded)

    for name, seq in msa_a2:
        df.loc[name, "name"] = name
        for pos in range(0, len(seq) - 100, 100):
            allele_finded = isIn(msa_a1, seq[pos: pos+400])
            if allele_finded:
                print(name, pos, allele_finded)
            df.loc[name, pos] = not (not allele_finded)

    df.to_csv(output_name, index=False)
    print(df)
    return df


def plot400bpMatchnessByPair(a1: str, a2: str) -> list[go.Figure]:
    """ deprecated: plot calculate400bpReadMatchness """
    index = "index/kir_2100_merge.save"
    df = pd.read_csv(f"{index}.400bp_match.{a1}_{a2}.csv")
    # df = df.fillna(-1)
    df = df.fillna(0)
    df_a1 = df[df.name.str.startswith(a1)]
    df_a2 = df[df.name.str.startswith(a2)]
    figs = []
    color = px.colors.qualitative.T10

    # a1
    data = np.array(df_a1.fillna(-1).iloc[:, 1:], dtype=int)
    figs.append(px.imshow(data,
                          title=f"{a1} mapped on {a2}",
                          aspect="auto",
                          color_continuous_scale=[color[-1], color[0], color[2]]))
    figs[-1].update_layout(xaxis={'ticktext': df.columns[1::5],
                                  'tickmode': "array",
                                  'tickvals': list(range(len(df.columns[1:])))[::5]})

    # a2
    data = np.array(df_a2.fillna(-1).iloc[:, 1:], dtype=int)
    figs.append(px.imshow(data,
                          title=f"{a2} mapped on {a1}",
                          aspect="auto",
                          color_continuous_scale=[color[-1], color[0], color[2]]))
    figs[-1].update_layout(xaxis={'ticktext': df.columns[1::5],
                                  'tickmode': "array",
                                  'tickvals': list(range(len(df.columns[1:])))[::5]})
    return figs


def plot400bpMatchness() -> list[go.Figure]:
    """ deprecated: Call calculate400bpReadMatchness with selected gene pairs """
    pairs = [
        ("KIR2DL1", "KIR2DS1"),
        ("KIR2DL1", "KIR2DL2"),
        ("KIR2DL1", "KIR2DL3"),
        ("KIR2DL5A", "KIR2DL5B"),
    ]
    for a1, a2 in pairs:
        calculate400bpReadMatchness(a1, a2)

    figs = []
    for a1, a2 in pairs:
        figs.extend(plot400bpMatchnessByPair(a1, a2))
    return figs


def plotVariantPositions() -> list[go.Figure]:
    """ Call plotVariantPosBetweenPairs for selected pairs """
    pairs = [
        ("KIR2DL1", "KIR2DS1"),
        ("KIR2DL2", "KIR2DL3"),
        ("KIR2DS3", "KIR2DS5"),
        ("KIR2DL5A", "KIR2DL5B"),
    ]
    figs = []
    for a1, a2 in pairs:
        figs.extend(plotVariantBetweenPair(a1, a2))
    return figs


def plotVariantBetweenPair(a1: str, a2: str) -> list[go.Figure]:
    """ plot variations between two gene along base """
    # read MSA
    index = "index/kir_2100_merge.save"
    msa = msaio.load_msa(f"{index}.KIR.fa", f"{index}.KIR.json")
    submsa = msa.select_allele(f"({a1})|({a2})").shrink()
    # calculate variantions
    base0 = submsa.get_variantion_base()
    base1 = submsa.select_allele(f"({a1})").get_variantion_base()
    base2 = submsa.select_allele(f"({a2})").get_variantion_base()
    base_between = set(base0) - set(base1) - set(base2)

    # plot
    figs = []
    title = f"Variantion of {a1} and {a2}"
    df = pd.DataFrame([
        *[{'pos': i, 'from': a1 + "'s variants"} for i in base1],
        *[{'pos': i, 'from': a2 + "'s variants"} for i in base2],
        *[{'pos': i, 'from': f"between {a1} and {a2}"} for i in base_between],
    ])
    figs.append(px.histogram(df, x='pos', color="from", nbins=100, title=title))
    return figs


def getVariantPosition(msa: Genemsa, allele_regex: str) -> set[int]:
    """ Find the variantions in the allele_regex selected alleles """
    return set(msa.select_allele(allele_regex).get_variantion_base())


def selectShortSequence(msa: Genemsa, minimal_threshold: float = 0.7) -> list[str]:
    allele_length = {name: len(seq.replace("-", "")) for name, seq in msa.alleles.items()}
    threshold = np.max(list(allele_length.values())) * minimal_threshold
    removed_alleles = []
    for allele, length in allele_length.items():
        if length < threshold:
            removed_alleles.append(allele)
    return removed_alleles


def calcVariantionsPosition() -> dict[tuple[str, str], set[int]]:
    """ calculate variantion position between geneA and geneB """
    # load data
    index = "index/kir_2100_merge.save"
    msa = msaio.load_msa(f"{index}.KIR.fa", f"{index}.KIR.json")
    output_name = f"{index}.gene_variation1"

    # read from saved data
    if os.path.exists(output_name + ".json"):
        diff_pos_json = json.load(open(output_name + ".json"))
        return {tuple(k.split(";")): set(v) for k, v in diff_pos_json.items()}  # type: ignore

    # Remove some strange sequences
    del msa.alleles['KIR*BACKBONE']

    # get All available gene
    genes = sorted(set(map(getGeneName, msa.alleles.keys())))

    # remove short sequences
    for gene in genes:
        msa_gene = msa.select_allele(f"{gene}\*")
        short_alleles = selectShortSequence(msa_gene)
        msa.remove(short_alleles)
        print("Remove", short_alleles)
        print(f"{gene} delete {len(short_alleles)} alleles in total {len(msa_gene.alleles)} allele")

    # get variantion across gene
    exes = {}
    diff_pos = {}
    with ProcessPoolExecutor(max_workers=threads) as executor:
        # run concurrent
        for gene1 in genes:
            for gene2 in genes:
                if gene1 <= gene2:
                    exes[(gene1, gene2)] = executor.submit(getVariantPosition,
                                                           msa, f"({gene1})|({gene2})")

        # gene1 vs gene2
        for key, exe in exes.items():
            diff_pos[key] = exe.result()

    # save
    json.dump({i + ";" + j: list(pos) for (i, j), pos in diff_pos.items()},
              open(output_name + ".json", "w"))
    return diff_pos


def plotDissimilarity() -> list[go.Figure]:
    """ plot dis-similarity matrix of all gene """
    diff_pos = calcVariantionsPosition()
    genes = set(chain.from_iterable(diff_pos.keys()))
    gene_id = {gene: id for id, gene in enumerate(sorted(genes))}

    # calculation
    diff_genes = []
    for gene1 in genes:
        for gene2 in genes:
            if gene1 == gene2:
                diff_genes.append({
                    'from': gene1,
                    'to': gene2,
                    'value': len(diff_pos[(gene1, gene2)]),
                })
            else:
                diff_genes.append({
                    'from': gene1,
                    'to': gene2,
                    'value': len(diff_pos[(min(gene1, gene2), max(gene1, gene2))]
                             - diff_pos[(gene1, gene1)]),
                })
    diff_genes_df = pd.DataFrame(diff_genes)
    print(diff_genes_df)
    print(diff_genes_df[diff_genes_df['value'] < 800])

    # pandas to 2d array
    data = np.zeros((len(genes), len(genes)))
    for d in diff_genes_df.to_dict(orient="records"):
        data[gene_id[d['from']], gene_id[d['to']]] = d['value']

    fig = px.imshow(data, text_auto=True, color_continuous_scale='RdBu_r',
                    width=1200, height=1200,
                    labels=dict(color="Base"),
                    x=list(gene_id.keys()), y=list(gene_id.keys())
                    )
    fig.update_layout(xaxis_side="top",
                      xaxis_title="from",
                      yaxis_title="to")

    return [fig]


def addGroupName(name: str) -> str:
    """
    Add addition read-group for each read by the gene generated from.

    This bam can be plotted on IGV and can be sorted by group name.

    Input: {name}.bam
    Return: {name}.rg.bam
    """
    outputname = name + ".rg"
    proc = runDocker("samtools", f"samtools view -h {name}.bam", capture_output=True)

    header = []
    data = []
    read_groups = set()
    for i in proc.stdout.split("\n"):
        if not i:
            continue
        if i.startswith("@"):
            header.append(i)
            continue
        rg = i.split('*')[0]
        data.append(i + f"\tRG:Z:{rg}")
        read_groups.add(rg)

    for rg in read_groups:
        header.append(f"@RG\tID:{rg}")

    with open(f"{outputname}.sam", "w") as f:
        f.writelines(map(lambda i: i + "\n", header))
        f.writelines(map(lambda i: i + "\n", data))

    samtobam(outputname)
    print(f"Save to {outputname}.bam")
    return ".rg"


def seqSplit(seq: list[str], sep: str) -> list[str]:
    """
    Split the seq with sep string

    Example:
        ```
        seq = ["abaa", "aababb"]
        sep = "ba"
        output: ["a", "baa", "aa", "babb"]
        ```
    """
    seqs = []
    for s in seq:
        s_split = s.split(sep)
        s_split[1:] = [sep + i for i in s_split[1:]]
        seqs.extend(s_split)
    return seqs


def selectKeyAlleleByResolution(allele_names: list[str],
                                resolution: str = "gene") -> list[str]:
    """
    Select alleles and all selected allels is different under the resolution

    Example:
        allele_names:
          - KIR2DL1*0010101
          - KIR2DL1*00102
          - KIR2DL1*00103
          - KIR2DL1*00201
          - KIR2DL1*0020201
          - KIR2DL1*003
        resolution: 3

        Return:
          - KIR2DL1*0010101
          - KIR2DL1*00201
          - KIR2DL1*003
    """
    s = set()
    reserved_name = []
    for i in sorted(allele_names):
        if resolution == "gene":
            key = i.split("*")[0]
        elif resolution == "7":
            key = i.split("*")[0] + "*" + i.split("*")[1][:7]
        elif resolution == "5":
            key = i.split("*")[0] + "*" + i.split("*")[1][:5]
        elif resolution == "5":
            key = i.split("*")[0] + "*" + i.split("*")[1][:3]
        else:
            key = i
        if key not in s:
            s.add(key)
            reserved_name.append(i)

    return sorted(reserved_name)


def extractRepeatSequence(index: str, gene: str) -> tuple[Genemsa, Genemsa]:
    """ Extract the KIR repeat region (about 500bp - 3000bp) """
    msa = msaio.load_msa(f"{index}.{gene}.fa",
                         f"{index}.{gene}.json")
    # get repeat seqs
    seq_list = [msa.get(f"{gene}*BACKBONE")[:5000]]
    seq_list = seqSplit(seq_list, "AATATGG")
    seq_list = seqSplit(seq_list, "ATAGGG")
    seq_list = seqSplit(seq_list, "GATA")
    seq_list = seqSplit(seq_list, "GATC")
    seq_list = seqSplit(seq_list, "GAGAGGAA")
    seq_list = seqSplit(seq_list, "GAGTGAG")
    seq_list = seqSplit(seq_list, "GATGTG")
    seq_list = seqSplit(seq_list, "GGTAG")
    seq_list = seqSplit(seq_list, "GGGAGCG")
    seq_list = seqSplit(seq_list, "GGGAGTGCG")
    seq_list = seqSplit(seq_list, "GGGGAGA")
    seq_list = seqSplit(seq_list, "GTAAT")
    seq_list = seqSplit(seq_list, "GTGAT")
    seq_list = seqSplit(seq_list, "GTTAT")

    # add position information along with sequence
    seq_pos_list = []
    pos = 0
    for seq in seq_list:
        seq_pos_list.append({
            'pos': pos,
            'seq': seq,
        })
        pos += len(seq)

    # reserve similar length (the pattern is from 17 to 22 bp)
    seq_pos_list = [i for i in seq_pos_list if 17 <= len(i['seq']) <= 22]  # type: ignore
    pprint(seq_pos_list)

    # run muscle
    file_repeat = f"{index}.{gene}.repeat"
    SeqIO.write([SeqRecord(Seq(i['seq']), id=f"{gene}*BB-{i['pos']:05d}")
                 for i in seq_pos_list],
                file_repeat + ".fa",
                "fasta")
    file_repeat = muscle(file_repeat)

    # print aligned repeat
    msa_repeat_element = Genemsa.from_MultipleSeqAlignment(
        AlignIO.read(file_repeat + ".fa", "fasta"))
    msa_repeat_element.append(f"{gene}*BB-con",
                              msa_repeat_element.get_consensus(include_gap=True))
    msa_repeat_element.set_reference(f"{gene}*BB-con")
    print(msa_repeat_element.sort_name().format_alignment_diff())
    return msa, msa_repeat_element


def displayRepeatRegions(msa: Genemsa,
                         msa_repeat_element: Genemsa,
                         resolution: str = "7") -> None:
    """ Display the MSA from extractRepeatSequence """
    # get repeat elements' start position and length
    repeat_regions = sorted([
        # pos                      length
        (int(name.split('-')[-1]), len(allele.replace("-", "")))
        for name, allele in msa_repeat_element.alleles.items()
        if "-con" not in name
    ])
    gene = msa_repeat_element.get_reference()[0].split("*")[0]

    # print possible sequences at each repeat regions
    for pos, length in repeat_regions:
        msa_part = msa[pos:pos+length]
        # aggregate same sequences
        seqs_count = Counter(msa_part.alleles.values())
        msa_part.alleles = {f"cluster{i:02d}-size-{c}": seq
                            for i, (seq, c) in enumerate(seqs_count.items())}
        print(msa_part.format_alignment_diff())

    # remove allele that are same under specific resolution
    alleles_to_show = selectKeyAlleleByResolution(msa.get_sequence_names(), resolution)
    msa = msa.select_allele(alleles_to_show)

    # print all sequences at each repeat regions
    for pos, length in repeat_regions:
        msa_part = msa[pos:pos+length]
        msa_part = msa_part.select_allele(alleles_to_show)
        print(msa_part.format_alignment_diff())


def plotAllelesLength(index: str) -> list[go.Figure]:
    """ plot all the length of allese in the dataset """
    kir = readFromMSAs(index)
    lengths = []
    for gene, msa in sorted(kir.items()):
        for allele, seq in msa.alleles.items():
            if "BACKBONE" in allele:
                continue
            lengths.append({
                'gene': gene,
                'allele': allele,
                'length': len(seq.replace('-', "")),
            })
    lengths_df = pd.DataFrame(lengths)

    summary = lengths_df.groupby("gene")["length"].describe()
    summary['mean'] = summary['mean'].apply(lambda i: "{:,.0f}".format(i))
    summary['std'] = summary['std'].apply(lambda i: "{:,.0f}".format(i))
    print(summary['mean'] + " ± " + summary['std'])
    return [px.box(lengths_df, x="gene", y="length", points="all")]


def statsAlleleNum(index: str) -> list[go.Figure]:
    """ plot all the length of allese in the dataset """
    kir = readFromMSAs(index)
    alleles = []
    for gene, msa in sorted(kir.items()):
        for allele, seq in msa.alleles.items():
            if "BACKBONE" in allele:
                continue
            alleles.append({
                'gene': gene,
                'allele': allele,
                'resolution': len(getAlleleField(allele, resolution=7)),
            })

    def resolutionCount(df: pd.DataFrame) -> pd.Series:
        count = {
            '>=3': sum(df["resolution"] >= 3),
            '>=5': sum(df["resolution"] >= 5),
            '>=7': sum(df["resolution"] >= 7),
        }
        return pd.Series(count)

    allele_df = pd.DataFrame(alleles)
    summary = allele_df.groupby("gene").apply(resolutionCount)
    print(summary)
    summary_plot = summary.stack().reset_index()
    summary_plot.columns = ["gene", "resolution", "count"]
    return [px.bar(summary_plot, x="gene", y="count", color="resolution", barmode="group")]


def calcVariationWithoutGap(msa: Genemsa) -> list[int]:
    """ Actually same as get_variantion_base but does not consider gap """
    freqs = msa.calculate_frequency()
    base = []
    for i, freq in enumerate(freqs):
        if sum(freq[:4]) not in freq[:4]:
            base.append(i)
    return base


def calcExonFlankingPosition(msa: Genemsa) -> list[list[int]]:
    """ Get all exon + flanking start and length (If overlap, merge it) """
    flanking_length = 200
    exon_position = []
    start_position = 0
    for block in msa.blocks:
        if block.type == "exon":
            exon_position.append((start_position, start_position + block.length))
        start_position += block.length

    # add flanking
    exon_with_flanking_position: list[list[int]] = []
    for a, b in exon_position:
        if (len(exon_with_flanking_position) and
                a - flanking_length <= exon_with_flanking_position[-1][1]):
            exon_with_flanking_position[-1][1] = min(b + flanking_length,
                                                     msa.get_length())
        else:
            exon_with_flanking_position.append([
                max(a - flanking_length, 0),
                min(b + flanking_length, msa.get_length())
            ])
    return exon_with_flanking_position


def addPercentage(df: pd.DataFrame) -> pd.DataFrame:
    """ Calculate the percentage of the value per gene """
    new_df = df.groupby(["gene", "region"]).first()
    percentage = new_df / df.groupby("gene").sum()
    new_df['percentage'] = percentage
    return new_df.reset_index()


def plotNumOfVariants() -> list[go.Figure]:
    """ Plot number of variant_num in exon and length of exon """
    # index = "index/kir_2100_merge.save"
    # msa = msaio.load_msa(f"{index}.KIR.fa", f"{index}.KIR.json")
    # genes = sorted(set(map(getGeneName, msa.alleles.keys())) - set(["KIR"]))
    # kir = {}
    # for gene in genes:
    #     kir[gene] = msa.select_allele(f"({gene})").shrink()
    index = "index/kir_2100_2dl1s1.save"
    index = "index/kir_2100_ab.save"
    index = "index/kir_2100_raw.save"
    kir = readFromMSAs(index)

    gene_variation_with_flanking_num = []
    gene_variation_num = []
    gene_length = []
    gene_length_with_flanking = []
    for gene, submsa in sorted(kir.items()):
        print(f"Process {gene}")
        # remove short alleles
        # short_alleles = selectShortSequence(submsa, 0.7)
        # submsa.remove(short_alleles)
        # print("Remove", short_alleles)
        # print(f"{gene} delete {len(short_alleles)} alleles in total {len(submsa.alleles)} allele")

        # calculate statistic
        gene_length.append({
            'gene': gene,
            'region': "exon",
            'length': submsa.select_exon().get_length(),
        })
        gene_length.append({
            'gene': gene,
            'region': "intron",
            'length': submsa.get_length() - gene_length[-1]['length'],  # type: ignore
        })

        exon_with_flanking_position = calcExonFlankingPosition(submsa)
        gene_length_with_flanking.append({
            'gene': gene,
            'region': "exon + 200bp flanking",
            'length': sum(map(lambda i: i[1] - i[0], exon_with_flanking_position)),
        })
        gene_length_with_flanking.append({
            'gene': gene,
            'region': "other",
            'length': submsa.get_length() - gene_length_with_flanking[-1]['length'],  # type: ignore
        })

        gene_variation_num.append({
            'gene': gene,
            'region': "exon",
            'variant_num': len(calcVariationWithoutGap(submsa.select_exon())),
        })
        gene_variation_num.append({
            'gene': gene,
            'region': "intron",
            'variant_num': len(calcVariationWithoutGap(submsa))
                           - gene_variation_num[-1]['variant_num'],  # type: ignore
        })

        msa_exon_with_flanking = reduce(lambda i, j: i + j, [submsa[a:b] for a, b in exon_with_flanking_position])
        gene_variation_with_flanking_num.append({
            'gene': gene,
            'region': "exon + 200bp flanking",
            'variant_num': len(calcVariationWithoutGap(msa_exon_with_flanking)),
        })
        gene_variation_with_flanking_num.append({
            'gene': gene,
            'region': "others",
            'variant_num': len(calcVariationWithoutGap(submsa))
                           - gene_variation_with_flanking_num[-1]['variant_num'],  # type: ignore
        })

    df_variant_num               = pd.DataFrame(gene_variation_num)
    df_variant_with_flanking_num = pd.DataFrame(gene_variation_with_flanking_num)
    df_length                    = pd.DataFrame(gene_length)
    df_length_with_flanking      = pd.DataFrame(gene_length_with_flanking)

    df_variant_num               = addPercentage(df_variant_num)
    df_variant_with_flanking_num = addPercentage(df_variant_with_flanking_num)
    df_length                    = addPercentage(df_length)
    df_length_with_flanking      = addPercentage(df_length_with_flanking)
    print(df_variant_num)
    print(df_variant_with_flanking_num)
    print(df_length)
    print(df_length_with_flanking)

    return [
        px.bar(df_variant_num,               x="gene", y="variant_num", color="region",
               text_auto="d",   title="Number of variants (exclude deletion)"),
        px.bar(df_variant_num,               x="gene", y="percentage",  color="region",
               text_auto=".0%", title="Number of variants (exclude deletion)"),
        px.bar(df_length,                    x="gene", y="length",      color="region",
               text_auto="d",   title="Length of Sequence"),
        px.bar(df_length,                    x="gene", y="percentage",  color="region",
               text_auto=".0%", title="Length of Sequence"),
        px.bar(df_variant_with_flanking_num, x="gene", y="variant_num", color="region",
               text_auto="d",   title="Number of variants (exclude deletion)"),
        px.bar(df_variant_with_flanking_num, x="gene", y="percentage",  color="region",
               text_auto=".0%", title="Number of variants (exclude deletion)"),
        px.bar(df_length_with_flanking,      x="gene", y="length",      color="region",
               text_auto="d",   title="Length of Sequence"),
        px.bar(df_length_with_flanking,      x="gene", y="percentage",  color="region",
               text_auto=".0%", title="Length of Sequence"),
    ]


if __name__ == "__main__":
    # displayRepeatRegions(*extractRepeatSequence("index/kir_2100_raw.save", "KIR3DL3"), "7")
    # displayRepeatRegions(*extractRepeatSequence("index/kir_2100_merge.save", "KIR"), "5")
    # addGroupName("index/kir_2100_merge.save.KIR")
    # addGroupName("index/kir_2100_merge_assign1.save.KIR")
    # exit()
    figs: list[go.Figure] = []
    # figs.extend(plot400bpMatchness())  # deprecated
    # figs.extend(plotVariantPositions())
    # figs.extend(plotDissimilarity())
    # figs.extend(plotNumOfVariants())
    # figs.extend(plotAllelesLength("index/kir_2100_ab.save"))
    figs.extend(statsAlleleNum("index5/kir_2100_split"))
    figs.extend(statsAlleleNum("index5/kir_2100_withexon_split"))

    # dash
    app = Dash(__name__)
    app.layout = html.Div([dcc.Graph(figure=f) for f in figs])
    app.run_server(debug=True, port=8052)
