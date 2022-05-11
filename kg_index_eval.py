import re
import numpy as np
import pandas as pd
import plotly.express as px
from dash import Dash, dcc, html, Input, Output
from pyHLAMSA import Genemsa
from concurrent.futures import ProcessPoolExecutor, as_completed


def calculate_400bp_read_matchness(a1, a2):
    # read
    msa = Genemsa.load_msa(f"index/kir_2100_merge.KIR.save.fa", "index/kir_2100_merge.KIR.save.json")

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


    df.to_csv(f"read_400bp_match_{a1}_{a2}.csv", index=False)
    print(df)


def plot_400bp_matchness(a1, a2):
    df = pd.read_csv(f"read_400bp_match_{a1}_{a2}.csv")
    # df = df.fillna(-1)
    df = df.fillna(0)
    df_a1 = df[df.name.str.startswith(a1)]
    df_a2 = df[df.name.str.startswith(a2)]
    figs = []
    color = px.colors.qualitative.T10

    # a1
    data = np.array(df_a1.fillna(-1).iloc[:, 1:], dtype=int)
    figs.append(px.imshow(data, title=a1, text_auto=True, aspect="auto", color_continuous_scale=[color[-1], color[0], color[2]]))
    figs[-1].update_layout(xaxis={'ticktext': df.columns[1::5], 'tickmode': "array", 'tickvals': list(range(len(df.columns[1:])))[::5]})

    # a2
    data = np.array(df_a2.fillna(-1).iloc[:, 1:], dtype=int)
    figs.append(px.imshow(data, title=a2, text_auto=True, aspect="auto", color_continuous_scale=[color[-1], color[0], color[2]]))
    figs[-1].update_layout(xaxis={'ticktext': df.columns[1::5], 'tickmode': "array", 'tickvals': list(range(len(df.columns[1:])))[::5]})

    return figs


def find_400bp_matchness():
    a1, a2 = "KIR2DL1", "KIR2DS1"
    a1, a2 = "KIR2DL5A", "KIR2DL5B"
    a1, a2 = "KIR2DL1", "KIR2DL2"
    calculate_400bp_read_matchness(a1, a2)
    return plot_400bp_matchness(a1, a2)


def diffBetweenAllele(msa, title):
    import plotly.express as px
    figs = []
    return figs[0]


def plot_variant_matchness():
    figs = []
    msa = Genemsa.load_msa(f"index/kir_2100_merge.KIR.save.fa", f"index/kir_2100_merge.KIR.save.json")

    gene = "KIR2DS1"
    base_gene = "KIR2DL1"
    # gene = "KIR2DL2"
    # base_gene = "KIR2DL3"
    # gene = "KIR2DS3"
    # base_gene = "KIR2DS5"

    # All genes
    # for g in set(map(lambda i: i.split("*")[0], msa.alleles.keys())):
    for g in [gene, base_gene]:
        submsa = msa.select_allele(f"({g})|({base_gene})").shrink()
        bs = submsa.get_variantion_base()
        print(g)
        print("Total length", submsa.get_length())
        print("Total base diff", len(bs))
        df = pd.DataFrame(bs, columns = ['pos'])
        figs.append( px.histogram(df, x='pos', title=f"{base_gene} vs {g}", width=600, height=800) )
        # figs.append( px.histogram(df, x='pos', title=title, width=600, height=800, histnorm='probability').update_layout(yaxis_tickformat = '.2%') )

    # diff variant of two = all variant - variant1 - variant2
    submsa = msa.select_allele(f"({gene})|({base_gene})").shrink()
    submsa_base0 = submsa.get_variantion_base()
    submsa_base1 = submsa.select_allele(f"({gene})").get_variantion_base()
    submsa_base2 = submsa.select_allele(f"({base_gene})").get_variantion_base()
    bases = set(submsa_base0) - set(submsa_base1) - set(submsa_base2)

    # plot diff base
    print(submsa.format_variantion_base())

    # plot
    title = f"Different between {base_gene} and {gene} (exclude internal variant)"
    print(title, len(bases))
    # .update_layout(yaxis_range=[0,70])
    df = pd.DataFrame(bases, columns = ['pos'])
    figs.append(px.histogram(df, x='pos', title=title, width=600, height=800))
    # px.histogram(df, x='pos', title=title, width=600, height=800, histnorm='probability').update_layout(yaxis_tickformat = '.2%') ])

    return figs


def getVariantPosition(msa, reg):
    return set(msa.select_allele(reg).get_variantion_base())


def plot_variant_across_gene():
    # load data
    msa = Genemsa.load_msa("index/kir_2100_merge.KIR.save.fa",
                           "index/kir_2100_merge.KIR.save.json")
    del msa.alleles['KIR*BACKBONE']
    del msa.alleles['KIR2DS4*0010103']
 
    # remove 2DL5 short sequences
    for i in map(lambda i: i[0], filter(lambda i: len(i[1].replace("-", "")) < 3000, msa.select_allele("KIR2DL5A.*").alleles.items())):
        print(f"delete {i}")
        del msa.alleles[i]
    for i in map(lambda i: i[0], filter(lambda i: len(i[1].replace("-", "")) < 3000, msa.select_allele("KIR2DL5B.*").alleles.items())):
        print(f"delete {i}")
        del msa.alleles[i]
 
    # get All available gene
    genes = sorted(set(map(lambda i: i.split("*")[0], msa.alleles.keys())))
    # genes.remove("KIR3DL3")
    # genes.remove("KIR3DP1")
    # genes.remove("KIR2DL5A")
    # genes.remove("KIR2DL5B")
    gene_id = {gene: id for id, gene in enumerate(genes)}

    # get variantion in gene
    bases = {}
    for gene in genes:
        bases[gene] = set(msa.select_allele(f"{gene}").get_variantion_base())

    # get variantion across gene
    diff_genes = []
    exes = {}
    with ProcessPoolExecutor(max_workers=8) as executor:
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

    # pandas to 2d array
    data = np.zeros((len(genes), len(genes)))
    for d in diff_genes:
        data[gene_id[d['from']], gene_id[d['to']]] = d['value']

    fig = px.imshow(data, text_auto=True, color_continuous_scale='RdBu_r',
                    width=1200, height=1200,
                    labels=dict(color="Base"),
                    x=list(gene_id.keys()), y=list(gene_id.keys())
          ).update_xaxes(side="top")

    # dash cannot show text on image
    # fig.write_image("genes_diff_all.png")
    # fig.write_image("genes_diff.png")  # no 3DP1
    # fig.write_image("genes_diff_noshort.png")  # no 3DP1 + no 2DL5 short
    fig.write_image("results/variant_diff_across_gene_noshort.png")  # no 3DP1 + no 2DL5 short
    return [fig]


if __name__ == "__main__":
    figs = []
    figs.extend(find_400bp_matchness())
    # figs.extend(plot_variant_matchness())
    # figs.extend(plot_variant_across_gene())

    # dash
    app = Dash(__name__)
    app.layout = html.Div([dcc.Graph(figure=f) for f in figs])
    app.run_server(debug=True, port=8052)
