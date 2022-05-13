import os
import numpy as np
import pandas as pd
import plotly.express as px
from dash import Dash, dcc, html, Input, Output
from pyHLAMSA import Genemsa
from concurrent.futures import ProcessPoolExecutor, as_completed
from kg_utils import threads, runDocker, runShell, samtobam


def calculate400bpReadMatchness(a1, a2):
    # read
    index = "index/kir_2100_merge.save"
    output_name = f"{index}.400bp_match.{a1}_{a2}.csv"
    if os.path.exists(output_name):
        return pd.read_csv(output_name)
    msa = Genemsa.load_msa(f"{index}.KIR.fa", f"{index}.KIR.json")

    # separate
    msa = msa.select_allele(f"{a1}.*|{a2}.*").shrink().reset_index().sort_name()
    msa_a1 = [(i[0], i[1].replace("-", "")) for i in msa.select_allele(f"{a1}.*").alleles.items()]
    msa_a2 = [(i[0], i[1].replace("-", "")) for i in msa.select_allele(f"{a2}.*").alleles.items()]

    def isIn(items, seq):
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


def plot400bpMatchness(a1, a2):
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
    figs.append(px.imshow(data, title=f"{a1} mapped on {a2}", aspect="auto", color_continuous_scale=[color[-1], color[0], color[2]]))
    figs[-1].update_layout(xaxis={'ticktext': df.columns[1::5], 'tickmode': "array", 'tickvals': list(range(len(df.columns[1:])))[::5]})

    # a2
    data = np.array(df_a2.fillna(-1).iloc[:, 1:], dtype=int)
    figs.append(px.imshow(data, title=f"{a2} mapped on {a1}", aspect="auto", color_continuous_scale=[color[-1], color[0], color[2]]))
    figs[-1].update_layout(xaxis={'ticktext': df.columns[1::5], 'tickmode': "array", 'tickvals': list(range(len(df.columns[1:])))[::5]})

    return figs


def plot400bpMatchness():
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
        figs.extend(plot400bpMatchness(a1, a2))
    return figs


def diffBetweenAllele(msa, title):
    # TODO what is this?
    pass


def plotVariantPosBetweenPairs():
    figs = []
    figs.extend(plotVariantPosBetweenPair("KIR2DS1", "KIR2DL1"))
    figs.extend(plotVariantPosBetweenPair("KIR2DL2", "KIR2DL3"))
    figs.extend(plotVariantPosBetweenPair("KIR2DS3", "KIR2DS5"))
    figs.extend(plotVariantPosBetweenPair("KIR2DL5A", "KIR2DL5"))
    return figs


def plotVariantPosBetweenPair(a1, a2):
    figs = []
    index = "index/kir_2100_merge.save"
    msa = Genemsa.load_msa(f"{index}.KIR.fa", f"{index}.KIR.json")

    """
    # gene vs base difference
    # for g in set(map(lambda i: i.split("*")[0], msa.alleles.keys())):
    for g in [gene, base_gene]:
        submsa = msa.select_allele(f"({g})|({base_gene})").shrink()
        title = f"SNP distribution of {base_gene} and {g}"
        bs = submsa.get_variantion_base()
        print(g)
        print("Total length", submsa.get_length())
        print("Total base diff", len(bs))
        df = pd.DataFrame(bs, columns = ['pos'])
        figs.append( px.histogram(df, x='pos', title=title, width=600, height=800) )
        # figs.append( px.histogram(df, x='pos', title=title, width=600, height=800, histnorm='probability').update_layout(yaxis_tickformat = '.2%') )
    """

    # Inter gene difference
    # diff variant of two = all variant - variant1 - variant2
    submsa = msa.select_allele(f"({a1})|({a2})").shrink()
    submsa_base0 = submsa.get_variantion_base()
    submsa_base1 = submsa.select_allele(f"({a1})").get_variantion_base()
    submsa_base2 = submsa.select_allele(f"({a2})").get_variantion_base()
    bases = set(submsa_base0) - set(submsa_base1) - set(submsa_base2)

    title = f"Internal Variantion inside {a1} and {a2}"
    df = pd.DataFrame([
        *[{'pos': i, 'from': a1} for i in submsa_base1],
        *[{'pos': i, 'from': a2} for i in submsa_base2],
    ])
    figs.append(px.histogram(df, x='pos', color="from", nbins=100, title=title, width=600, height=800))

    # plot
    title = f"Different between {a1} and {a2} (exclude internal variant)"
    print(title, len(bases))
    # .update_layout(yaxis_range=[0,70])
    df = pd.DataFrame(bases, columns=['pos'])
    figs.append(px.histogram(df, x='pos', title=title, nbins=100, width=600, height=800))
    # px.histogram(df, x='pos', title=title, width=600, height=800, histnorm='probability').update_layout(yaxis_tickformat = '.2%') ])

    return figs


def getVariantPosition(msa, reg):
    return set(msa.select_allele(reg).get_variantion_base())


def plotVariationAllGene():
    # load data
    index = "index/kir_2100_merge.save"
    msa = Genemsa.load_msa(f"{index}.KIR.fa", f"{index}.KIR.json")

    # Remove some strange sequences
    del msa.alleles['KIR*BACKBONE']
    del msa.alleles['KIR2DS4*0010103']
    # remove 2DL5 short sequences
    for i in map(lambda i: i[0], filter(lambda i: len(i[1].replace("-", "")) < 3000, msa.select_allele("KIR2DL5.*").alleles.items())):
        print(f"delete {i}")
        del msa.alleles[i]

    # get All available gene
    genes = sorted(set(map(lambda i: i.split("*")[0], msa.alleles.keys())))
    # genes.remove("KIR3DL3")
    # genes.remove("KIR3DP1")
    # genes.remove("KIR2DL5A")
    # genes.remove("KIR2DL5B")
    gene_id = {gene: id for id, gene in enumerate(genes)}

    outputfile = f"{index}.variation.csv"
    if not os.path.exists(outputfile):
        # get variantion in gene
        bases = {}
        for gene in genes:
            bases[gene] = set(msa.select_allele(f"{gene}").get_variantion_base())

        # get variantion across gene
        diff_genes = []
        exes = {}
        with ProcessPoolExecutor(max_workers=threads) as executor:
            # run concurrent
            for gene1 in genes:
                for gene2 in genes:
                    if gene1 != gene2:
                        exes[(gene1, gene2)] = executor.submit(getVariantPosition, msa, f"({gene1})|({gene2})")

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
        diff_genes_df.to_csv(outputfile, index=False)
    else:
        diff_genes_df = pd.read_csv(outputfile)

    # pandas to 2d array
    data = np.zeros((len(genes), len(genes)))
    for d in diff_genes_df.to_dict(orient="records"):
        data[gene_id[d['from']], gene_id[d['to']]] = d['value']

    fig = px.imshow(data, text_auto=True, color_continuous_scale='RdBu_r',
                    width=1200, height=1200,
                    labels=dict(color="Base"),
                    x=list(gene_id.keys()), y=list(gene_id.keys())
                    ).update_xaxes(side="top")

    return [fig]


def addGroupName(name):
    """ Input {name}.bam Return {name}.rg.bam """
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
        else:
            rg = i.split('*')[0]
            data.append(i + f"\tRG:Z:{rg}")
            read_groups.add(rg)

    for rg in read_groups:
        header.append(f"@RG\tID:{rg}")

    with open(f"{outputname}.sam", "w") as f:
        f.writelines(map(lambda i: i + "\n", header))
        f.writelines(map(lambda i: i + "\n", data))

    samtobam(outputname)
    print("Save to {outputname}.bam")
    return ".rg"


if __name__ == "__main__":
    # addGroupName("index/kir_2100_merge.save.KIR")
    figs = []
    # figs.extend(plot400bpMatchness())
    # figs.extend(plotVariantPosBetweenPairs())
    figs.extend(plotVariationAllGene())

    # dash
    app = Dash(__name__)
    app.layout = html.Div([dcc.Graph(figure=f) for f in figs])
    app.run_server(debug=True, port=8052)
