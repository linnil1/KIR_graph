import os
import re
import json
from pprint import pprint
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor
from scipy.stats import norm
import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# plot
from dash import Dash, dcc, html
import plotly.graph_objects as go
import plotly.express as px
# import dash_bio as dashbio

from kg_utils import getSamples
from kg_eval import EvaluateKIR
from kg_typping import readLocus, readVariants, isInExon

typing_by = "linnil1"
count_only_exon = False
# linnil1 hisat


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def collectAlleleName(reads_alleles):
    """ Get all alleles from all the reads """
    allele_names = set()
    for read_alleles in reads_alleles:
        for typ, variants_alleles in read_alleles.items():
            if typ not in ["lp", "ln", "rp", "rn"]:
                continue
            for variant_alleles in variants_alleles:
                allele_names.update(set(variant_alleles))
    return allele_names


def flatten(its):
    for it in its:
        yield from it


def variant2onehot(allele_name_map, variant_alleles):
    """ Make possible allele into onehot encoding inside a variant """
    a = np.zeros(len(allele_name_map), dtype=bool)
    for allele in variant_alleles:
        a[allele_name_map[allele]] = True
    return a


def onehot2Prob(onehot, pos=True):
    """ onehot encoding array to prob """
    # TODO: use quality
    prob = np.ones(onehot.shape)
    if pos:
        prob[onehot] = 0.999
        prob[np.logical_not(onehot)] = 0.001
    else:
        prob[onehot] = 0.001
        prob[np.logical_not(onehot)] = 0.999
    return prob


def getReadProb(allele_name_map, read_alleles):
    """
    Calculate the probility of alleles that can generate the read

    Input size: variants_in_reads x alleles_num
    Return size: 1 x alleles_num
    """
    if getNH(read_alleles['l_sam']) != 1:
        return None

    def variant2prob(variant_alleles, pos=True):
        return onehot2Prob(variant2onehot(allele_name_map, variant_alleles), pos=pos)
    prob = [
        *[variant2prob(variant_alleles, True)  for variant_alleles in read_alleles['lp']],
        *[variant2prob(variant_alleles, True)  for variant_alleles in read_alleles['rp']],
        *[variant2prob(variant_alleles, False) for variant_alleles in read_alleles['ln']],
        *[variant2prob(variant_alleles, False) for variant_alleles in read_alleles['rn']],
    ]
    if not prob:  # strnage case?
        return None
    prob = np.stack(prob).prod(axis=0)

    # too ambiguous
    # if not any(prob >= 1e-6):
    #     return None

    if test_alleles:
        p_allele = set(flatten(flatten([read_alleles['lp'], read_alleles['rp']])))
        n_allele = set(flatten(flatten([read_alleles['ln'], read_alleles['rn']])))
        test_id = np.array(list(allele_name_map[i] for i in test_alleles))
        read_id = read_alleles["l_sam"].split("\t")[0]
        # print allele but called negative
        # if "2DL1" not in read_id and "2DS1" not in read_id:
        #     return prob
        #     return None
        read_allele = read_id.split('-')[0]
        if not allele_name_map.get(read_allele):
            return prob
        if prob[allele_name_map[read_allele]] <= 1e-3:
            ok = True
            for alleles, v in zip(read_alleles['ln'], read_alleles['lnv']):
                if read_allele in alleles:
                    print(alleles, v)
                    ok = False
            if not ok:
                print("left", read_alleles['l_sam'])
            ok = True
            for alleles, v in zip(read_alleles['rn'], read_alleles['rnv']):
                if read_allele in alleles:
                    print(alleles, v)
                    ok = False
            if not ok:
                print("right", read_alleles['r_sam'])
        """
        # print the loss compare to other
        if prob[test_id[:len(test_alleles)//2]].max() != prob[test_id[len(test_alleles)//2:]].max():
        # if prob[test_id[:len(test_alleles)//2]].max() < 1e-3:
            print(prob[test_id], "\n",
                 'pos', [allele for allele in p_allele if allele in test_alleles], "\n",
                 'neg', [allele for allele in n_allele if allele in test_alleles],
                 read_alleles,
            )
        """
    return prob


def printProb(data):
    """
2
0 -52.37587690948907
  id   1 name KIR3DS1*0130108      fraction 0.35294117647058826
  id   2 name KIR3DS1*078          fraction 0.6470588235294118
1 -52.37587690948907
  id   0 name KIR3DS1*0130107      fraction 0.4215686274509804
  id   2 name KIR3DS1*078          fraction 0.5784313725490197
    """
    n = data['n']
    for idx, candidate in enumerate(zip(data['allele_id'],
                                        data['allele_name'],
                                        data['fraction'],
                                        data['fraction_norm'],
                                        data['unique_count'],
                                        data['value'])):
        if idx >= 10:
            break
        print(idx, candidate[-1])
        for allele_id, allele_name, frac, frac_norm, unique_count in zip(*candidate[:-1]):
            print(f"  id {allele_id:3} name {allele_name:20s}"
                  f" fraction {frac:.5f} "
                  f" unique_count {unique_count}")


def calcProb(log_probs, prob_last, prob_iter=1, top_n=300, allele_length={}):
    """
    Calculate the probility of
    allele_1 x allele_2 x ... x allele_n
    can generated the reads with max probility
    """
    assert prob_iter >= 1
    if prob_iter == 1:
        # find the maximum probility of allele across all reads
        # remove_n = int(len(log_probs) * 0.01)
        # remove_n = 20
        # prob_1 = np.sort(log_probs, axis=0)[remove_n:].sum(axis=0)
        prob_1 = log_probs.sum(axis=0)
        prob_1_index = np.array(list(reversed(np.argsort(prob_1)[-top_n:])))
        # additional value
        prob_1_top = log_probs[:, prob_1_index]
        prob_1_top_allele = np.array([[i] for i in prob_1_index])
        # N = reads
        return {
            'n': 1,
            'value': prob_1[prob_1_index],                                                      # top_n
            'allele_id': prob_1_top_allele,                                                     # top_n x n
            'best_prob': prob_1_top,                                                            # top_n x n
            'fraction': np.ones(prob_1_top_allele.shape),                                       # top_n x n
            'fraction_norm': np.ones(prob_1_top_allele.shape) / allele_length[prob_1_top_allele],  # top_n x n
            'unique_count': np.ones(prob_1_top_allele.shape),                                   # top_n x n
        }

    def uniqueAllele(alleles, indexs):
        """
        Input: m * allele_n
        Output: after_unique_m * allele_n
        """
        unique_allele = []
        unique_indexs = []
        s = set()
        for ind in indexs:
            allele = alleles[ind]
            allele = tuple(sorted(allele))
            if allele in s:
                continue
            s.add(allele)
            unique_allele.append(list(allele))
            unique_indexs.append(ind)
        return np.array(unique_allele), np.array(unique_indexs)

    # get previous top-n allele
    prob_1_top = prob_last['best_prob']
    prob_1_top_allele = prob_last['allele_id']

    # Find the maximum in
    # (top-n x (allele-1 ... allele-n-1)) x (allele_1 ... allele_m)
    prob_2 = np.maximum(log_probs, prob_1_top.T[:, :, None])
    prob_2 = prob_2.sum(axis=1).flatten()
    # prob_2 = np.sort(prob_2, axis=1)[:, remove_n:].sum(axis=1).flatten()
    # Find the mean
    # prob_2 = (log_probs + prob_1_top.T[:, :, None] * (prob_iter - 1)).sum(axis=1).flatten() / prob_iter
    prob_2_allele = np.hstack([
        # [[1],[3],[5]] -> [[1],[3],[5],[1],[3],[5]]
        np.repeat(prob_1_top_allele, log_probs.shape[1], axis=0),
        # [1,3,5] -> [1,3,5,1,3,5]
        np.tile(np.arange(log_probs.shape[1]), len(prob_1_top_allele))[:, None],
    ])
    prob_2_index = np.array(list(reversed(np.argsort(prob_2)[-top_n:])))

    # additional value
    prob_2_top_allele, prob_2_index = uniqueAllele(prob_2_allele, prob_2_index)
    prob_2_top = log_probs[:, prob_2_top_allele].max(axis=2)
    prob_2_value = prob_2[prob_2_index]
    prob_2_belong = np.equal(log_probs[:, prob_2_top_allele], prob_2_top[:, :, None])  # reads x top-n x allele_num(0/1)

    # Count of uniquely assigned
    prob_2_unique_mapped = prob_2_belong.sum(axis=2) == 1
    prob_2_unique_count = prob_2_belong.copy()
    prob_2_unique_count[np.logical_not(prob_2_unique_mapped)] = 0
    prob_2_unique_count = prob_2_unique_count.sum(axis=0)

    prob_2_fraction = (prob_2_belong / prob_2_belong.sum(axis=2)[:, :, None]).sum(axis=0)
    prob_2_fraction = prob_2_fraction / prob_2_fraction.sum(axis=1, keepdims=True)
    prob_2_fraction_norm = prob_2_fraction / allele_length[prob_2_top_allele]
    # disable
    # prob_2_fraction_nrom = prob_2_fraction

    # sort by loss + fraction (evenly -> better)
    value_fraction = np.vstack([-prob_2_value, np.abs(prob_2_fraction_norm - prob_2_fraction_norm.mean(axis=1, keepdims=True)).sum(axis=1)]).T.tolist()
    rank_value_fraction = sorted(range(len(value_fraction)), key=lambda i: value_fraction[i])
    prob_2_top = prob_2_top[:, rank_value_fraction]
    prob_2_top_allele = prob_2_top_allele[rank_value_fraction]
    prob_2_value = prob_2_value[rank_value_fraction]
    prob_2_fraction = prob_2_fraction[rank_value_fraction]
    prob_2_fraction_norm = prob_2_fraction_norm[rank_value_fraction]
    prob_2_unique_count = prob_2_unique_count[rank_value_fraction]

    return {
        'n': prob_iter,
        'allele_id': prob_2_top_allele,
        'best_prob': prob_2_top,
        'value': prob_2_value,
        'fraction': prob_2_fraction,
        'fraction_norm': prob_2_fraction_norm,
        'unique_count': prob_2_unique_count,
    }


def getReadProbTmp(allele_name_map, reads_alleles_chunk):
    return [getReadProb(allele_name_map, read_alleles) for read_alleles in reads_alleles_chunk]


def calcProbInit(allele_names, reads_alleles):
    allele_length_map = readAlleleLength(index)
    # name -> id
    allele_name_map = dict(zip(allele_names, range(len(allele_names))))
    allele_name_map_rev = {b: a for a, b in allele_name_map.items()}
    allele_length = np.array([allele_length_map[allele_name_map_rev[i]] for i in range(len(allele_names))])
    allele_length = allele_length / 10000.

    # prob per read

    if typing_by != "merge":
        probs = getReadProbTmp(allele_name_map, reads_alleles)
    else:
        print("Alleles in read -> vector")
        with ProcessPoolExecutor(max_workers=30) as executor:
            exes = []
            probs = []
            for i in range(0, len(reads_alleles), 400):
                exes.append(executor.submit(getReadProbTmp, allele_name_map, reads_alleles[i:i+400]))
            for exe in tqdm(exes):
                prob = exe.result()
                probs.extend(prob)
        print("Alleles in read -> vector Done")

    # reads x alleles
    probs = [i for i in probs if i is not None]
    probs = np.stack(probs)

    norm_probs = probs / probs.sum(axis=1, keepdims=True)
    log_probs = np.log10(probs)
    # log_probs = np.log10(norm_probs)

    # debug
    if test_alleles:
        test_id = np.array(list(allele_name_map[i] for i in test_alleles))
        print(log_probs[:, test_id[:len(test_alleles)//2]].max(axis=1).sum())
        print(log_probs[:, test_id[len(test_alleles)//2:]].max(axis=1).sum())

    """
    # plot prob
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
    """
    return log_probs, allele_length, allele_name_map_rev


def calcProbAll(allele_names, reads_alleles, iter_max=5):
    """
    Calculate the probility of iter_max of alleles
    to fit the positive and negative variants in each reads
    """
    # init
    log_probs, allele_length, allele_name_map_rev = calcProbInit(allele_names, reads_alleles)
    top_n = 300  # reduce computation
    prob_save = []

    # find the maximum probility of allele across all reads
    for prob_iter in range(1, 1 + iter_max):
        if prob_save:
            prob_last = prob_save[-1]
        else:
            prob_last = None
        prob_res = calcProb(log_probs, prob_last,
                            allele_length=allele_length,
                            prob_iter=prob_iter,
                            top_n=top_n)
        prob_res['allele_name'] = [list(map(allele_name_map_rev.get, i)) for i in prob_res['allele_id']]
        prob_save.append(prob_res)
    return prob_save


def readAlleleLength(index):
    """ Calculate the length of all alleles """
    seq_len = {}
    for seq in SeqIO.parse(f"{index}_sequences.fa", "fasta"):
        seq_len[seq.id] = len(seq.seq)
    return seq_len


def avgGeneLength(seq_len):
    """ Average the length of all alleles inside a gene """
    gene_length = defaultdict(list)
    for k, v in seq_len.items():
        gene_length[getGeneName(k)].append(v)
    gene_length = {k: np.mean(v) for k, v in gene_length.items()}
    # merge 2DL5A 2DL5B
    gene_length["KIR2DL5"] = (gene_length["KIR2DL5A"] + gene_length["KIR2DL5B"]) * 0.5
    # merge 2DL1 2DS1
    gene_length["KIR2DL1S1"] = (gene_length["KIR2DL1"] + gene_length["KIR2DS1"]) * 0.5
    return gene_length


def getGeneName(s):
    return s.split("*")[0]


def getAlleleName(s):
    return s.split("-")[0]


def getNH(sam_info):
    """ Extract NH from record. Note: It return 1 / NH """
    NH_re = re.findall(r"NH:i:(\d+)", sam_info)
    if NH_re:
        return int(NH_re[0])
    else:
        return 1


def plotAlignConfusion(data, level="gene", weighted=True, no_multi=False):
    allele_length = readAlleleLength("index/kir_2100_raw.mut01")
    gene_length = avgGeneLength(allele_length)

    # gene count
    gene_from_to_count = defaultdict(int)
    for gene, reads_alleles in data.items():
        for read_alleles in reads_alleles:
            record = read_alleles['l_sam']
            if no_multi and getNH(record) != 1:
                continue
            if level == "gene":
                id = (getGeneName(record), getGeneName(gene))
            elif level == "allele":
                id = (getAlleleName(record), getGeneName(gene))
            if weighted:
                gene_from_to_count[id] += 1 / getNH(record)
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
    for ref_name, gene_data in report_data.items():
        print(ref_name)
        for i, c in enumerate(gene_data['count'][:10]):
            print("  ", i, *c)
        for i, p in enumerate(gene_data['prob']):
            if p[1] < 0.001:
                continue
            print("  ", i, *p)
    return report_data


def typingWithReportAndCopyNumber(report, gene_cn):
    called_alleles = []
    for ref_name, cn in gene_cn.items():
        if not cn:
            continue
        est_prob = 1 / cn
        alleles = report[ref_name]['prob'][:cn]
        # called_alleles.extend([i[0] for i in alleles])
        # TODO:
        # case1: a1=0.66 a2=0.33
        # case2: a1=0.9 a2=0.1
        for allele_name, allele_prob in alleles:
            pred_count = max(1, round(allele_prob / est_prob))
            for i in range(min(cn, pred_count)):
                called_alleles.append(allele_name)
            cn -= pred_count
            if cn <= 0:
                break

    return called_alleles


def hasAllele(variants_alleles, allele_name):
    """ Is str in list[list[str]] """
    for variant_alleles in variants_alleles:
        if allele_name in variant_alleles:
            return True
    return False


def plotPosNegRate(reads_alleles, predict_alleles):
    count = defaultdict(lambda: {'all': 0, 'pos': 0, 'neg': 0, 'both': 0, 'non': 0})
    for read_alleles, predict_allele in zip(reads_alleles, predict_alleles):
        # print(read_alleles)
        count[predict_allele]['all'] += 1
        pos = hasAllele(read_alleles['rp'], predict_allele) or hasAllele(read_alleles['lp'], predict_allele)
        neg = hasAllele(read_alleles['rn'], predict_allele) or hasAllele(read_alleles['ln'], predict_allele)
        if pos and neg:
            count[predict_allele]['both'] += 1
        elif pos and not neg:
            count[predict_allele]['pos'] += 1
        elif not pos and neg:
            count[predict_allele]['neg'] += 1
        else:
            count[predict_allele]['non'] += 1
    pprint(count)
    return count


def plotPosNegRateWithReadAns(data):
    counts = {}
    for gene, reads_alleles in data.items():
        ans_alleles = [getAlleleName(read_alleles['l_sam']) for read_alleles in reads_alleles]
        print(gene, "ans", len(ans_alleles))
        print(gene, "all", len(reads_alleles))
        count = plotPosNegRate(reads_alleles, ans_alleles)
        for allele_name, d in count.items():
            if allele_name.startswith(getGeneName(gene)):
                counts[allele_name] = d
    pprint(counts)


def plotConfustion(df, title=""):
    figs = []
    order = {'from': sorted(set(df['from'])), 'to': sorted(set(df['to']))}
    figs.append(px.bar(df, x="from", y="value",           color="to",   text='to',   category_orders=order, title=title))
    figs.append(px.bar(df, x="to",   y="value",           color="from", text='from', category_orders=order,))
    figs.append(px.bar(df, x="from", y="norm_value_from", color="to",   text='to',   category_orders=order,))
    figs.append(px.bar(df, x="to",   y="norm_value_to",   color="from", text='from', category_orders=order,))
    return figs


norm_max = 2.0


def getCNDist(base, base_dev=0.02):
    """
    Return ( CN x norm_read_count(500bin) ) array
    indicate the probility of read_count belong to the CN
    """
    x = np.linspace(0, norm_max, 500)
    space = norm_max / 500
    cn = np.arange(1, 10)
    # y0 = norm.pdf(x, loc=base*0.1, scale=base_dev*5)
    y0 = norm.pdf(x, loc=0, scale=base_dev*1.5)
    y = np.stack([y0, *[norm.pdf(x, loc=base*(n), scale=base_dev*n) for n in cn]])
    return x, y * space


def getCNDistBest(norm_read_count, assume_base=None):
    loss = []  # loss
    num_read_count, _ = np.histogram(norm_read_count, bins=500, range=(0, norm_max))
    space = norm_max / 500
    base_dev = 0.02

    for base in np.arange(0, norm_max, 0.01):
        _, y = getCNDist(base, base_dev)
        max_cn_prob = y.max(axis=0)  # We will choose the highest probility of CN at that read count
        loss.append((base, np.sum( np.log(max_cn_prob + 1e-9) * num_read_count )))

    loss = np.array(loss)
    max_point = loss[np.argmax(loss[:, 1]), :]
    base = max_point[0]
    if assume_base and base / assume_base > 1.7:  # TODO: maybe more CN?
        base /= 2
    x, y = getCNDist(base, base_dev)
    max_cn_arg = y.argmax(axis=0)
    return {
        'base': base,
        'x': x,
        'dist': y,
        'loss': loss,
        'count': norm_read_count,
        'class': [max_cn_arg[int(cn / space)] for cn in norm_read_count]
    }


def plotCNDist(result):
    fig_loss = px.line(x=result['loss'][:, 0], y=result['loss'][:, 1])
    fig_loss.add_vline(x=result['base'], line_dash="dash", line_color="gray",
                       annotation_text=f"max value={result['base']}")
    fig_loss.update_layout(
        xaxis_title="read-depth / 300",
        yaxis_title="likelihood",
    )

    fig_dist = go.Figure()
    fig_dist.add_vline(x=result['base'], line_dash="dash", line_color="gray",
                       annotation_text=f"CN=1 value={result['base']}")
    fig_dist.add_trace(go.Scatter(x=result['count'],
                                  y=[0 for i in result['count']],
                                  mode='markers', name="Observed"))
    fig_dist.update_layout(
        xaxis_title="read-depth / 300",
        yaxis_title="probility",
    )
    for n in range(len(result['dist'])):
        x = result['x']
        if n < 6:
            fig_dist.add_trace(go.Scatter(x=x, y=result['dist'][n], name=f"cn={n}"))
    return [fig_loss, fig_dist]


def getCNPerRead(data):
    allele_length = readAlleleLength(index)
    gene_length = avgGeneLength(allele_length)

    gene_read_count = defaultdict(int)
    for gene, reads_alleles in data.items():
        for read_alleles in reads_alleles:
            gene_read_count[gene] += 1 / getNH(read_alleles['l_sam'])
        gene_read_count[gene] = gene_read_count[gene] / gene_length[getGeneName(gene)]
    gene_read_count = list(gene_read_count.items())
    print("Normalized read count")
    pprint(gene_read_count)

    norm_read_count = [i[1] for i in gene_read_count]
    result = getCNDistBest(norm_read_count)


    return (
        dict(zip([i[0] for i in gene_read_count], result['class'])),
        plotCNDist(result),
    )


def typingWithCopyNumber(data, gene_cn):
    called_alleles = []
    for gene, reads_alleles in data.items():
        # if gene != "KIR2DS1*BACKBONE":
        #     continue
        allele_names = sorted(collectAlleleName(reads_alleles))
        # remove multiple aligned
        # this written in getReadProb
        # reads_alleles = [read_allele for read_allele in reads_alleles if getNH(read_allele['l_sam']) == 1]
        if not gene_cn.get(gene):
            continue

        prob_save = calcProbAll(allele_names, reads_alleles, iter_max=gene_cn[gene])
        if test_alleles:
            for prob in prob_save:
                print(prob['n'])
                printProb(prob)

        data = prob_save[gene_cn[gene] - 1]

        # TODO: testing the new selection
        best_id = 0
        for i in range(len(data['fraction'])):
            expect_prob = 1 / gene_cn[gene]
            if all(j >= expect_prob / 2 for j in data['fraction'][i]):
                break
            best_id += 1
        if best_id >= len(data['fraction']):
            best_id = 0

        called_alleles.extend(data['allele_name'][best_id])
        print(f"{gene}: CN={gene_cn[gene]}, alleles={data['allele_name'][0]} score={data['value'][0]}")
        # print(data['fraction'][0])
    return called_alleles


def getCNbyVariant(index, data, assume_depth=None):
    # insertion deletion has more count of course
    snp = pd.read_csv(f"{index}.snp", sep="\t", header=None)
    single_id = set(snp.loc[snp[1] == "single", 0])

    variants = []
    if count_only_exon:
        # ugly method
        exon_region = readLocus(index)
        snp_variants = readVariants(index)
        for v in snp_variants:
            v.in_exon = isInExon(exon_region[v.ref], v)
        snp_variants = {v.id: v for v in snp_variants}

    for gene, reads_alleles in data.items():
        # if gene != "KIR3DP1*BACKBONE":
        #     continue
        pos = []
        neg = []
        for read_alleles in reads_alleles:
            if getNH(read_alleles['l_sam']) != 1:
                continue
            for key in ['lpv', 'rpv']:
                pos.extend(read_alleles[key])
            for key in ['lnv', 'rnv']:
                neg.extend(read_alleles[key])
        pos = filter(lambda i: snp_variants[i].in_exon, pos)
        neg = filter(lambda i: snp_variants[i].in_exon, neg)
        pos = Counter(pos)
        neg = Counter(neg)

        for v in (pos.keys() | neg.keys()) & single_id:
            pos_num = pos.get(v, 0)
            neg_num = neg.get(v, 0)
            num = pos_num + neg_num
            # if num <= 20:
            #     continue
            variants.append({
                'ratio': pos_num / num,
                'pos_num': pos_num,
                'neg_num': neg_num,
                'total': num,
                'name': v,
                'gene': gene,
            })

    df = pd.DataFrame(variants)
    figs = []
    order_total = [i['name'] for i in sorted(variants, key=lambda i: i['total'])]
    figs.append( px.bar(df, x="name", y=["pos_num", "neg_num"],
                        category_orders={'name': order_total},
                        title="Depth of variants") )
    figs.append( px.bar(df, x="name", y=["total"], color="gene",
                        category_orders={'name': order_total},
                        color_discrete_sequence=px.colors.qualitative.Dark24,
                        title="Depth of variants group by gene") )
    order_ratio = [i['name'] for i in sorted(variants, key=lambda i: i['ratio'])]
    figs.append( px.bar(df[df['ratio'] > 0.05], x="name", y="ratio", color="gene",
                        color_discrete_sequence=px.colors.qualitative.Dark24,
                        category_orders={'name': order_ratio},
                        title="Ratio of variants") )
    figs.append( px.box(df, x="gene", y="total",
                        category_orders={'gene': sorted(data.keys())},
                        title="Depth of each gene (no del/ins)") )
    if not count_only_exon:
        gene_total = df.groupby(by="gene", as_index=False)['total'].median()
    else:
        gene_total = df.groupby(by="gene", as_index=False)['total'].quantile(.75)
    print(gene_total)

    # copy from before
    # 300 is a magic number (to Scale)
    norm_read_count = list(gene_total['total'] / 300)
    if not assume_depth:
        result = getCNDistBest(norm_read_count)
    else:
        result = getCNDistBest(norm_read_count, assume_base=assume_depth / 300)
    figs.extend(plotCNDist(result))
    return dict(zip(gene_total['gene'], result['class'])), figs


def getCNbySamtoolsDepth(index, name, assume_depth=None):
    df = pd.read_csv(f"{name}.no_multi.depth.tsv", sep="\t", header=None, names=["gene", "pos", "depth"])
    figs = []
    figs.append(px.line(df, x="pos", y="depth", color="gene", title="Read Depth"))

    if count_only_exon:
        exon_region = readLocus(index)
        df_exon = []
        for gene, region in exon_region.items():
            for (start, end) in region:
                df_exon.append(df[ (df["gene"] == gene) & (start <= df["pos"]) & (df["pos"] <= end) ])
        df_exon = pd.concat(df_exon)
        df = df_exon

    # count
    figs.append(px.box(df, x="gene", y="depth", title="Read Depth"))
    if not count_only_exon:
        gene_total = df.groupby(by="gene", as_index=False)['depth'].median()
    else:
        gene_total = df.groupby(by="gene", as_index=False)['depth'].quantile(.75)
    print(gene_total)

    # same as before
    # 300 is a magic number (to Scale)
    norm_read_count = list(gene_total['depth'] / 300)
    if not assume_depth:
        result = getCNDistBest(norm_read_count)
    else:
        result = getCNDistBest(norm_read_count, assume_base=assume_depth / 300)
    figs.extend(plotCNDist(result))
    # pprint(dict(zip(gene_total['gene'], result['class'])))
    return dict(zip(gene_total['gene'], result['class'])), figs


def typingUntilEnd(index, data, assume_depth=None):
    links = {}
    with open(f"{index}.link") as f:
        for i in f:
            if not i:
                continue
            id, alleles = i.strip().split("\t")
            links[id] = alleles.split(' ')

    for gene, reads_alleles in data.items():
        assert gene == "KIR*BACKBONE"
        allele_names = sorted(set(flatten(links.values())))
        for read_alleles in reads_alleles:
            read_alleles['lp'] = [links[v] for v in read_alleles['lpv']]
            read_alleles['ln'] = [links[v] for v in read_alleles['lnv']]
            read_alleles['rp'] = [links[v] for v in read_alleles['rpv']]
            read_alleles['rn'] = [links[v] for v in read_alleles['rnv']]
        log_probs, allele_length, allele_name_map_rev = calcProbInit(allele_names, reads_alleles)
        top_n = 100  # reduce computation
        prob_save = []

        # find the maximum probility of allele across all reads
        for prob_iter in range(1, 27 + 1):
            if prob_save:
                prob_last = prob_save[-1]
            else:
                prob_last = None
            prob_res = calcProb(log_probs, prob_last,
                                allele_length=allele_length,
                                prob_iter=prob_iter,
                                top_n=top_n)
            prob_res['allele_name'] = [list(map(allele_name_map_rev.get, i)) for i in prob_res['allele_id']]
            prob_save.append(prob_res)

            print(prob_res['n'])
            printProb(prob_res)
        json.dump(prob_save, open("test_save_prob.json", "w"), cls=NumpyEncoder)
    return []


def typingPerSample(index, name, assume_depth=None):
    if typing_by == "hisat":
        hisat_report = readReport(name)
    print("Read json")
    data = json.load(open(name + ".json"))
    print("Read json done")
    figs = []
    # align_gene_confustion_df = plotAlignConfusion(data, no_multi=True)
    # figs.extend(plotConfustion(align_gene_confustion_df, title=f"Reads (gene level, no-duplication) {name}"))
    '''
    align_gene_confustion_df = plotAlignConfusion(data)
    figs.extend(plotConfustion(align_gene_confustion_df, title=f"Reads (gene level, weighted) {name}"))
    align_gene_confustion_df = plotAlignConfusion(data, weighted=False)
    figs.extend(plotConfustion(align_gene_confustion_df, title="Reads (gene level)"))
    align_gene_confustion_df = plotAlignConfusion(data, level="allele")
    figs.extend(plotConfustion(align_gene_confustion_df, title="Reads (allele level)"))
    # TODO: 2DL5
    '''
    if typing_by == "merge":
        return typingUntilEnd(index, data, assume_depth=assume_depth), figs
    # plotPosNegRateWithReadAns(data)
    # gene_cn, fig = getCNPerRead(data)
    gene_cn, fig = getCNbySamtoolsDepth(index, name, assume_depth=assume_depth)
    # gene_cn, fig = getCNbyVariant(index, data, assume_depth=assume_depth)
    figs.extend(fig)
    # gene_cn = {"KIR2DL1S1*BACKBONE": 3}
    pprint(gene_cn)
    if typing_by == "linnil1":
        called_alleles = typingWithCopyNumber(data, gene_cn)
    elif typing_by == "hisat":
        called_alleles = typingWithReportAndCopyNumber(hisat_report, gene_cn)
    else:
        raise NotImplementedError()
    # print(called_alleles)

    """
    if test_alleles:
        for gene, reads_alleles in data.items():
            for allele in test_alleles:
                plotPosNegRate(reads_alleles, [allele for read_alleles in reads_alleles]))
    """

    return called_alleles, figs


if __name__ == "__main__":
    depth = None
    # basename = "data/linnil1_syn_wide.test10"
    # index = "index/kir_2100_raw.mut01"
    # ans_csv = "linnil1_syn_wide/linnil1_syn_wide.summary.csv"
    basename = "data/linnil1_syn_30x"
    basename = "data/linnil1_syn_exon"
    index = "index/kir_2100_2dl1s1.mut01"
    suffix = ".index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr"
    # index = "index/kir_2100_merge.mut01"
    # suffix = ".index_kir_2100_merge.mut01.hisatgenotype.errcorr"
    # typing_by = "merge"
    ans_csv = ""
    if "linnil1_syn_30x" in basename:
        depth = 30
        ans_csv = "linnil1_syn_30x/linnil1_syn_30x.summary.csv"
    if "linnil1_syn_exon" in basename:
        ans_csv = "linnil1_syn_exon/linnil1_syn_exon.summary.csv"
        depth = 90
        count_only_exon = True
    if "linnil1_syn_wide" in basename:
        ans_csv = "linnil1_syn_wide/linnil1_syn_wide.summary.csv"

    names = getSamples(basename, suffix + ".json", return_name=True)
    names = names[0:10]
    figs = []
    called_alleles_dict = {}
    dfs = []

    test_alleles = []
    # test_alleles = ["KIR3DL3*0030104", "KIR3DL3*0090101", "KIR3DL3*0030104", "KIR3DL3*00801"]
    # test_alleles = ["KIR2DL1*0030230", "KIR2DS1*0020103", "KIR2DS1*0020115", "KIR2DL1*0030230", "KIR2DS1*0020104", "KIR2DS1*013"]
    # test_alleles = ["KIR2DS5*019", "KIR2DS5*019", "KIR2DS5*019", "KIR2DS5*0270102"]
    # test_alleles = ["KIR2DL1*0020112", "KIR2DL1*0030228", "KIR2DS1*0020108", "KIR2DS1*0020115"]
    # test_alleles += ["KIR2DL1*0020112", "KIR2DL1*0030228", "KIR2DS1*0020111", "KIR2DS1*013"]

    for name, name_id in names:
        name = os.path.splitext(name)[0]
        print(name, name_id)
        called_alleles, fig = typingPerSample(index, name, assume_depth=depth)
        figs.extend(fig)
        called_alleles_dict[name_id] = called_alleles
        df = pd.DataFrame([{
            'id':  name_id,
            'name': name,
            'alleles': "_".join(called_alleles),
        }])
        new_suffix = f".type_by_{typing_by}"
        # df.to_csv(name + f".linnil1{new_suffix}.csv", index=False, sep="\t")
        dfs.append(df)
    dfs = pd.concat(dfs, ignore_index=True)
    print(dfs)

    if ans_csv:
        EvaluateKIR(ans_csv).compareCohert(called_alleles_dict, skip_empty=True)
    app = Dash(__name__)
    app.layout = html.Div([dcc.Graph(figure=f) for f in figs])
    app.run_server(debug=True, port=8051)
