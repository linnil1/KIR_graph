import re
import json
from pprint import pprint
from collections import defaultdict
from scipy.stats import norm
import numpy as np
import pandas as pd
from Bio import SeqIO

# plot
from dash import Dash, dcc, html
import plotly.graph_objects as go
import plotly.express as px
import dash_bio as dashbio


def collectAlleleName(reads_alleles):
    allele_names = set()
    for read_alleles in reads_alleles:
        for typ, variants_alleles in read_alleles.items():
            if typ not in ["lp", "ln", "rp", "rl"]:
                continue
            for variant_alleles in variants_alleles:
                allele_names.update(set(variant_alleles))
    return allele_names


def variant2onehot(allele_name_map, variant_alleles):
    a = np.zeros(len(allele_name_map), dtype=int)
    for allele in variant_alleles:
        a[allele_name_map[allele]] = 1
    return a


def onehot2Prob(onehot, pos=True):
    prob = np.ones(onehot.shape)
    if pos:
        prob[onehot] = 0.999
        prob[np.logical_not(onehot)] = 0.001
    else:
        prob[onehot] = 0.001
        prob[np.logical_not(onehot)] = 0.999
    return prob


def getReadProb(allele_name_map, read_alleles):
    def variant2prob(variant_alleles, pos=True):
        return onehot2Prob(variant2onehot(allele_name_map, variant_alleles), pos=pos)
    prob = [
        *[variant2prob(variant_alleles, True)  for variant_alleles in read_alleles['lp']],
        *[variant2prob(variant_alleles, True)  for variant_alleles in read_alleles['rp']],
        *[variant2prob(variant_alleles, False) for variant_alleles in read_alleles['ln']],
        *[variant2prob(variant_alleles, False) for variant_alleles in read_alleles['rn']],
    ]
    if not prob:  # strnage case?
        return np.ones(len(allele_name_map))

    prob = np.stack(prob).prod(axis=0)
    return prob


def printProb(data):
    n = data['n']
    for idx, candidate in enumerate(zip(data['allele_id'],
                                        data['allele_name'],
                                        data['fraction'],
                                        data['value'])):
        if idx >= 20:
            break
        print(idx, candidate[-1])
        for allele_id, allele_name, frac in zip(*candidate[:-1]):
            print(f"  id {allele_id:3} name {allele_name:20s} fraction {frac}")


def uniqueAllele(alleles):
    unique_alleles = []
    s = set()
    for allele in alleles:
        allele = tuple(sorted(allele))
        if allele in s:
            continue
        s.add(allele)
        unique_alleles.append(list(allele))
    return np.array(unique_alleles)


def calcProb(allele_names, reads_alleles):
    allele_name_map = dict(zip(allele_names, range(len(allele_names))))
    allele_name_map_rev = {b: a for a, b in allele_name_map.items()}

    # prob per read
    probs = []
    for read_alleles in reads_alleles:
        prob = getReadProb(allele_name_map, read_alleles)
        if prob is not None:
            probs.append(prob)
    probs = np.stack(probs)
    norm_probs = probs / probs.sum(axis=1, keepdims=True)
    log_probs = np.log(norm_probs)
    # log_probs = np.log(probs)

    # init
    prob_save = []
    top_n = 30    # reduce computation
    prob_iter_max = 5

    prob_1 = log_probs.sum(axis=0)
    # print(prob_1.shape, prob_1)
    prob_1_index = np.array(list(reversed(np.argsort(prob_1)[-top_n:])))
    # print(prob_1_index.shape, prob_1_index)
    prob_1_top = log_probs[:, prob_1_index]
    # print(prob_1_top.shape, prob_1_top)
    prob_1_top_allele = np.array([[i] for i in prob_1_index])
    # print(prob_1_top_allele.shape, prob_1_top_allele)
    # print(prob_1.shape, prob_1_top.shape, prob_1_top_allele.shape)

    # N = reads
    prob_save.append({
        'n': 1,
        'allele_id': prob_1_top_allele,                                                     # top_n x n
        'allele_name': [list(map(allele_name_map_rev.get, i)) for i in prob_1_top_allele],  # top_n x n
        'best_prob': prob_1_top,                                                            # top_n x N
        'value': prob_1[prob_1_index],                                                      # top_n
        'fraction': np.ones(prob_1_top_allele.shape),                                       # top_n x n
    })

    for prob_iter in range(2, 1 + prob_iter_max):
        prob_1_top = prob_save[-1]['best_prob']
        prob_1_top_allele = prob_save[-1]['allele_id']

        prob_2 = np.maximum(log_probs, prob_1_top.T[:, :, None]).sum(axis=1).flatten()
        # print(prob_2.shape)
        prob_2_allele = np.hstack([
            np.repeat(prob_1_top_allele, log_probs.shape[1], axis=0),
            np.tile(np.arange(log_probs.shape[1]), len(prob_1_top_allele))[:, None],
        ])
        # print(prob_2_allele.shape)
        prob_2_index = np.array(list(reversed(np.argsort(prob_2)[-top_n:])))
        # print(prob_2_index.shape)
        prob_2_top_allele = uniqueAllele(prob_2_allele[prob_2_index])

        # print(prob_2_top_allele.shape)
        prob_2_top = log_probs[:, prob_2_top_allele].max(axis=2)
        # print(prob_2_top.shape)
        prob_2_value = prob_2[prob_2_index]
        # print(prob_2_value.shape)
        prob_2_belong = np.equal(log_probs[:, prob_2_top_allele], prob_2_top[:, :, None])
        # print(prob_2_belong.shape)
        prob_2_fraction = (prob_2_belong / prob_2_belong.sum(axis=2)[:, :, None]).sum(axis=0)
        prob_2_fraction = prob_2_fraction / prob_2_fraction.sum(axis=1, keepdims=True)
        # print(prob_2_fraction.shape)

        prob_save.append({
            'n': prob_iter,
            'allele_id': prob_2_top_allele,
            'allele_name': [list(map(allele_name_map_rev.get, i)) for i in prob_2_top_allele],
            'best_prob': prob_2_top,
            'value': prob_2_value,
            'fraction': prob_2_fraction,
        })
    return prob_save


def plotConfustion(df, title=""):
    figs = []
    order = {'from': sorted(set(df['from'])), 'to': sorted(set(df['to']))}
    figs.append(px.bar(df, x="from", y="value",           color="to",   text='to',   category_orders=order, title=title))
    figs.append(px.bar(df, x="to",   y="value",           color="from", text='from', category_orders=order,))
    figs.append(px.bar(df, x="from", y="norm_value_from", color="to",   text='to',   category_orders=order,))
    figs.append(px.bar(df, x="to",   y="norm_value_to",   color="from", text='from', category_orders=order,))
    return figs


def readAlleleLength(index):
    seq_len = {}
    for seq in SeqIO.parse(f"{index}_sequences.fa", "fasta"):
        seq_len[seq.id] = len(seq.seq)
    return seq_len


def avgGeneLength(seq_len):
    gene_length = defaultdict(list)
    for k, v in seq_len.items():
        gene_length[getGeneName(k)].append(v)
    gene_length = {k: np.mean(v) for k, v in gene_length.items()}
    gene_length["KIR2DL5"] = (gene_length["KIR2DL5A"] + gene_length["KIR2DL5B"]) * 0.5
    return gene_length


def getGeneName(s):
    return s.split("*")[0]


def getAlleleName(s):
    return s.split("-")[0]


def getNH(sam_info):
    NH_re = re.findall(r"NH:i:(\d+)", sam_info)
    if NH_re:
        return 1 / int(NH_re[0])
    else:
        return 1


def plotAlignConfusion(data, level="gene", weighted=True):
    allele_length = readAlleleLength("index/kir_2100_raw.mut01")
    gene_length = avgGeneLength(allele_length)

    # gene count
    gene_from_to_count = defaultdict(int)
    for gene, reads_alleles in data.items():
        for read_alleles in reads_alleles:
            record = read_alleles['l_sam']
            if level == "gene":
                id = (getGeneName(record), getGeneName(gene))
            elif level == "allele":
                id = (getAlleleName(record), getGeneName(gene))
            if weighted:
                gene_from_to_count[id] += getNH(record)
            else:
                gene_from_to_count[id] += 1

    # create dataframe
    gene_from_to = []
    for (gfrom, gto), count in gene_from_to_count.items():
        gene_from_to.append({'from': gfrom,
                             'to':   gto,
                             'from_length': gene_length[gfrom] if level == "gene" else allele_length[gfrom],
                             'to_length':   gene_length[gto],
                             'value': count,
                            })
    df = pd.DataFrame(gene_from_to)
    df['norm_value_from'] = df['value'] / df['from_length']
    df['norm_value_to']   = df['value'] / df['to_length']
    return df


def readReport(name):
    report_data = json.load(open(name + ".report.json"))
    for gene_data in report_data:
        count, prob = gene_data
        for i, c in enumerate(count[:10]):
            print(i, c)
        for i, p in enumerate(prob):
            if p[1] < 0.001:
                continue
            print(i, p)


def hasAllele(variants_alleles, allele_name):
    for variant_alleles in variants_alleles:
        if allele_name in variant_alleles:
            return True
    return False


def plotPosNegRate(data):
    count = defaultdict(lambda: {'all': 0, 'pos': 0, 'neg': 0, 'both': 0, 'non': 0})
    for gene, reads_alleles in data.items():
        for read_alleles in reads_alleles:
            # print(read_alleles)
            allele_name = getAlleleName(read_alleles['l_sam'])
            count[allele_name]['all'] += 1
            pos = hasAllele(read_alleles['rp'], allele_name) or hasAllele(read_alleles['lp'], allele_name)
            neg = hasAllele(read_alleles['rn'], allele_name) or hasAllele(read_alleles['ln'], allele_name)
            if pos and neg:
                count[allele_name]['both'] += 1
                print(read_alleles)
            elif pos and not neg:
                count[allele_name]['pos'] += 1
            elif not pos and neg:
                count[allele_name]['neg'] += 1
                print(read_alleles)
            else:
                count[allele_name]['non'] += 1
    pprint(count)


def getCNDist(base, base_dev=0.02):
    """
    Return ( CN x norm_read_count(500bin) ) array
    indicate the probility of read_count belong to the CN
    """
    x = np.linspace(0, 1.5, 500)
    space = 1.5 / 500
    cn = np.arange(1, 10)
    y0 = norm.pdf(x, loc=base*0.2, scale=base_dev*5)
    y = np.stack([y0, *[norm.pdf(x, loc=base*(n+0.2), scale=base_dev*n) for n in cn]])
    return y * space


def getCNDistBest(norm_read_count):
    loss = []  # loss
    num_read_count, _ = np.histogram(norm_read_count, bins=500, range=(0, 1.5))
    space = 1.5 / 500

    for base in np.arange(0, 1, 0.01):
        y = getCNDist(base)
        max_cn_prob = y.max(axis=0)  # We will choose the highest probility of CN at that read count
        loss.append((base, np.sum( np.log(max_cn_prob) * num_read_count )))

    loss = np.array(loss)
    max_point = loss[np.argmax(loss[:, 1]), :]
    base = max_point[0]
    y = getCNDist(base)
    max_cn_arg = y.argmax(axis=0)
    return {
        'base': base,
        'dist': y,
        'loss': loss,
        'class': [max_cn_arg[int(cn / space)] for cn in norm_read_count]
    }


def getCNPerRead(data):
    allele_length = readAlleleLength(index)
    gene_length = avgGeneLength(allele_length)

    gene_read_count = defaultdict(int)
    for gene, reads_alleles in data.items():
        for read_alleles in reads_alleles:
            gene_read_count[gene] += getNH(read_alleles['l_sam'])
        gene_read_count[gene] = gene_read_count[gene] / gene_length[getGeneName(gene)]
    gene_read_count = list(gene_read_count.items())
    print("Normalized read count")
    pprint(gene_read_count)

    norm_read_count = [i[1] for i in gene_read_count]
    result = getCNDistBest(norm_read_count)

    fig_loss = px.line(x=result['loss'][:, 0], y=result['loss'][:, 1])
    fig_dist = go.Figure()
    fig_dist.add_trace(go.Scatter(x=norm_read_count,
                                  y=[0 for i in norm_read_count],
                                  mode='markers', name="Observed"))
    for n in range(len(result['dist'])):
        x = np.linspace(0, 1.5, 500)
        if n  < 5:
            fig_dist.add_trace(go.Scatter(x=x, y=result['dist'][n], name=f"cn={n}"))

    return (
        list(zip([i[0] for i in gene_read_count], result['class'])),
        [fig_loss, fig_dist],
    )


def extractAnswerFromSummary(id):
    data = pd.read_csv("linnil1_syn_wide/summary.csv", sep="\t", header=None)
    data = list(data[2].str.split("_"))
    return sorted(data[id])


def evaluteAllele(ans_list, predict_list):
    predict_list = sorted(predict_list)
    match_pair = []

    for gene in set([getGeneName(i) for i in ans_list + predict_list]):
        match_pair.extend(
            evaluteGeneAllele([i for i in ans_list if i.startswith(gene)],
                              [i for i in predict_list if i.startswith(gene)]))
    match_pair.sort(key=lambda i: i[1] or i[2])
    for t, a, b in match_pair:
        if t == "Match":
            print(f"{a:18} OK {b:18}")
        elif t == "Mismatch":
            print(f"{a:18} <3 {b:18}")
        elif t == "FN":
            print(f"{a:18} <-")
        elif t == "FP":
            print(f"{'':18} -> {b:18}")
    return match_pair


def evaluteGeneAllele(a, b):
    match_pair = []
    i = 0
    for allele in list(b):
        if allele in a:
            a.remove(allele)
            b.remove(allele)
            match_pair.append(("Match", allele, allele))

    for allele_b in list(b):
        for allele_a in a:
            if allele_a[:allele_a.index("*") + 4] == allele_b[:allele_b.index("*") + 4]:
                a.remove(allele_a)
                b.remove(allele_b)
                match_pair.append(("Mismatch", allele_a, allele_b))
                break

    for allele in a:
        match_pair.append(("FN", allele, None))
    for allele in b:
        match_pair.append(("FP", None, allele))
    return match_pair


def evaluateResult(results):
    TP, FP, FN, match_gene, match_3, total, cn_error = 0, 0, 0, 0, 0, 0, 0
    fail_allele, fail_sample = 0, 0
    for result in results:
        if all(i[0] == 'FN' for i in result):
            fail_allele += len(result)
            fail_sample += 1
            continue
        total   += len([i[1] for i in result if i[1]])
        TP      += sum([i[0] == "Match" for i in result])
        FN      += sum([i[0] == "FN" for i in result])
        FP      += sum([i[0] == "FP" for i in result])
        match_3 += sum([i[0] == "Mismatch" for i in result])
        ans_gene = [i[1].split("*")[0] for i in result if i[0] == "FN" ]
        prd_gene = [i[2].split("*")[0] for i in result if i[0] == "FP" ]
        for g in prd_gene:
            if g in ans_gene:
                ans_gene.remove(g)
                match_gene += 1
            else:
                cn_error += 1
        cn_error += len(ans_gene)

    print(f"Cannot called", f"sample={fail_sample} alleles={fail_allele}")
    print(f"Total alleles {total}")
    print(f"  * {TP=}")
    print(f"  * Match_3_digits={match_3} ")
    print(f"  * {FN=} {FP=}")
    print(f"    * Match_gene={match_gene}")
    print(f"    * Copy_number_error={cn_error}")


def typingWithCopyNumber(data, gene_cn):
    called_alleles = []
    for gene, reads_alleles in data.items():
        allele_names = sorted(collectAlleleName(reads_alleles))
        prob_save = calcProb(allele_names, reads_alleles)

        if gene_cn[gene]:
            data = prob_save[gene_cn[gene] - 1]
            called_alleles.extend(data['allele_name'][0])
            print(f"{gene}: CN={gene_cn[gene]}, alleles={data['allele_name'][0]} score={data['value'][0]}")
            # print(data['fraction'][0])
    return called_alleles


def typingPerSample(name):
    # readReport(name)
    data = json.load(open(name + ".json"))
    # data = json.load(open("test2DL5.json"))
    figs = []
    align_gene_confustion_df = plotAlignConfusion(data)
    figs.extend(plotConfustion(align_gene_confustion_df, title="Reads (gene level, weighted)"))
    '''
    align_gene_confustion_df = plotAlignConfusion(data, weighted=False)
    figs.extend(plotConfustion(align_gene_confustion_df, title="Reads (gene level)"))
    align_gene_confustion_df = plotAlignConfusion(data, level="allele")
    figs.extend(plotConfustion(align_gene_confustion_df, title="Reads (allele level)"))
    # TODO: 2DL5
    plotPosNegRate(data)
    '''
    gene_cn, fig = getCNPerRead(data)
    gene_cn = dict(gene_cn)
    pprint(gene_cn)
    figs.extend(fig)

    called_alleles = typingWithCopyNumber(data, gene_cn)
    print(called_alleles)

    '''
    for gene, reads_alleles in data.items():
        for prob in prob_save:
            print(prob['n'])
            printProb(prob)
        figs.extend([
            px.imshow(np.log(probs), text_auto=True, aspect="auto", title=gene),
            px.imshow(norm_probs ** 0.1, text_auto=True, aspect="auto"),
            # px.imshow(np.log(norm_probs), text_auto=True, aspect="auto"),
        ])
        figs.extend([dashbio.Clustergram(
            data=probs,
            column_labels=allele_names,
            hidden_labels='row',
            width=2000,
            height=1000,
            cluster="row",
        )])
    '''
    return called_alleles, figs


if __name__ == "__main__":
    index = "index/kir_2100_raw.mut01"
    name = "data/linnil1_syn_wide.00.kir_2100_raw.mut01.hisatgenotype"
    name_id = 0

    names = [name]
    name_ids = [name_id]
    figs = []
    called_alleles_evals = []

    for name, name_id in zip(names, name_ids):
        called_alleles, fig = typingPerSample(name)
        figs.extend(fig)
        called_alleles_evals.append(
            evaluteAllele(extractAnswerFromSummary(name_id), called_alleles)
        )

    evaluateResult(called_alleles_evals)

    app = Dash(__name__)
    app.layout = html.Div([dcc.Graph(figure=f) for f in figs])
    app.run_server(debug=True, port=8051)
