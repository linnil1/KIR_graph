import os
import json
import numpy as np
import pandas as pd
from pprint import pprint
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import norm
from dash import Dash, dcc, html
from collections import defaultdict, Counter
from Bio import SeqIO
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema

from kg_typping import Hisat2, getNH
from kg_eval import EvaluateKIR


def readSamtoolsDepth(depth_filename):
    df = pd.read_csv(depth_filename, sep="\t",
                     header=None, names=["gene", "pos", "depth"])
    return df


def cropSamtoolsDepthInExon(df, exon_region):
    df_exon = []
    for gene, region in exon_region.items():
        for (start, end) in region:
            df_exon.append(df[ (df["gene"] == gene) & (start <= df["pos"]) & (df["pos"] <= end) ])
    df_exon = pd.concat(df_exon)
    return df_exon



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


def reconstructReadsVariants(reads_filename):
    data = json.load(open(reads_filename))
    save_reads, id_map_variant = data['reads'], data['variants']
    for reads in save_reads.values():
        for read in reads:
            read['lp'] = [id_map_variant[v] for v in read['lpv']]
            read['ln'] = [id_map_variant[v] for v in read['lnv']]
            read['rp'] = [id_map_variant[v] for v in read['rpv']]
            read['rn'] = [id_map_variant[v] for v in read['rnv']]
    return save_reads


def removeMultipleAlignmentRead(save_reads):
    return {
        gene: [read for read in reads if getNH(read['l_sam']) == 1]
        for gene, reads in save_reads.items()
    }


def readReference(index):
    return SeqIO.to_dict(SeqIO.parse(f"{index}_backbone.fa", "fasta"))


def readReport(report_json_filename):
    report_data = json.load(open(report_json_filename))
    for ref_name, gene_data in report_data.items():
        print(ref_name)
        for i, c in enumerate(gene_data['count'][:10]):
            print("  ", i, *c)
        for i, p in enumerate(gene_data['prob']):
            if p[1] < 0.001:
                continue
            print("  ", i, *p)
    return report_data


def typingWithReport(report, cn):
    est_prob = 1 / cn
    alleles = report['prob'][:cn]
    # called_alleles.extend([i[0] for i in alleles])
    # TODO:
    # case1: a1=0.66 a2=0.33
    # case2: a1=0.9 a2=0.1
    called_alleles = []
    for allele_name, allele_prob in alleles:
        pred_count = max(1, round(allele_prob / est_prob))
        for i in range(min(cn, pred_count)):
            called_alleles.append(allele_name)
        cn -= pred_count
        if cn <= 0:
            break

    return called_alleles


def merge2DL5CN(gene_cn):
    new_gene = "KIR2DL5"
    genes = ["KIR2DL5A", "KIR2DL5B"]
    for gene in genes:
        if new_gene not in gene_cn:
            gene_cn[new_gene] = 0
        if gene in gene_cn:
            gene_cn[new_gene] += gene_cn.pop(gene)


def merge2DL1S1CN(gene_cn):
    new_gene = "KIR2DL1S1"
    genes = ["KIR2DL1", "KIR2DS1"]
    for gene in genes:
        if new_gene not in gene_cn:
            gene_cn[new_gene] = 0
        gene_cn[new_gene] += gene_cn.pop(gene)


def readAlleleResult(filename):
    df = pd.read_csv(filename, sep="\t")
    records = df.to_dict('records')
    for record in records:
        record['alleles'] = record['alleles'].split("_")
    return records


def readCNResult(filename):
    df = pd.read_csv(filename, sep="\t")
    cn = {}
    for record in df.itertuples():
        cn[record.gene] = record.cn
    return cn


class AlleleTypping:
    def __init__(self, reads, allele_length_map):
        self.top_n = 300
        self.print_num = 10

        self.allele_name_map = {}
        self.allele_name_map_rev = {}
        self.log_probs = None
        self.reads2AlleleProb(reads)

        allele_length = np.array([allele_length_map[self.allele_name_map_rev[i]] for i in range(len(self.allele_name_map_rev))])
        self.allele_length = allele_length / 10000.  # just a non-important normalize
        self.result = []

    @staticmethod
    def collectAlleleNames(backbones_reads):
        """ Get all alleles from all the reads """
        allele_names = set()
        for reads in backbones_reads:
            for typ, variants in reads.items():
                if typ not in ["lp", "ln", "rp", "rn"]:
                    continue
                for variant in variants:
                    allele_names.update(set(variant['allele']))
        return allele_names

    def alleles2onehot(self, alleles):
        """ Make possible allele into onehot encoding inside a variant """
        a = np.zeros(len(self.allele_name_map), dtype=bool)
        for allele in alleles:
            a[self.allele_name_map[allele]] = True
        return a

    @staticmethod
    def onehot2Prob(onehot):
        """ onehot encoding array to prob """
        # TODO: use quality
        prob = np.ones(onehot.shape)
        prob[onehot] = 0.999
        prob[np.logical_not(onehot)] = 0.001
        return prob

    def calcProbInRead(self, read):
        return [
            *[self.onehot2Prob(               self.alleles2onehot(variant['allele']) ) for variant in read['lp']],
            *[self.onehot2Prob(               self.alleles2onehot(variant['allele']) ) for variant in read['rp']],
            *[self.onehot2Prob(np.logical_not(self.alleles2onehot(variant['allele']))) for variant in read['ln']],
            *[self.onehot2Prob(np.logical_not(self.alleles2onehot(variant['allele']))) for variant in read['rn']],
        ]

    def reads2AlleleProb(self, reads):
        allele_names = self.collectAlleleNames(reads)
        self.allele_name_map = dict(zip(allele_names, range(len(allele_names))))
        self.allele_name_map_rev = {b: a for a, b in self.allele_name_map.items()}

        probs = []
        for read in reads:
            prob = self.calcProbInRead(read)
            if not len(prob):
                continue
            probs.append(np.stack(prob).prod(axis=0))

        probs = [i for i in probs if i is not None]
        probs = np.stack(probs)
        log_probs = np.log10(probs)
        self.probs = probs
        self.log_probs = log_probs

    @staticmethod
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

    def typingByNum(self, n):
        assert len(self.result) == 0
        for _ in range(n):
            self.addCandidateNum()
            self.printProb(-1)
        return self.result[-1]

    def addCandidateNum(self):
        if not len(self.result):
            # find the maximum probility of allele across all reads
            # remove_n = int(len(log_probs) * 0.01)
            # remove_n = 20
            # prob_1 = np.sort(log_probs, axis=0)[remove_n:].sum(axis=0)
            prob_1 = self.log_probs.sum(axis=0)
            prob_1_index = np.array(list(reversed(np.argsort(prob_1)[-self.top_n:])))
            # additional value
            prob_1_top = self.log_probs[:, prob_1_index]
            prob_1_top_allele = np.array([[i] for i in prob_1_index])

            self.result.append({
                'n': 1,
                'value': prob_1[prob_1_index],                                                      # top_n
                'allele_id': prob_1_top_allele,                                                     # top_n x n
                'best_prob': prob_1_top,                                                            # top_n x n
                'fraction': np.ones(prob_1_top_allele.shape),                                       # top_n x n
                'fraction_norm': np.ones(prob_1_top_allele.shape) / self.allele_length[prob_1_top_allele],  # top_n x n
                'unique_count': np.ones(prob_1_top_allele.shape),                                   # top_n x n
                'allele_name': [list(map(self.allele_name_map_rev.get, i)) for i in prob_1_top_allele],     # top_n x n
            })
            return self.result[-1]

        # get previous top-n allele
        prob_1_top = self.result[-1]['best_prob']
        prob_1_top_allele = self.result[-1]['allele_id']

        # Find the maximum in
        # (top-n x (allele-1 ... allele-n-1)) x (allele_1 ... allele_m)
        prob_2 = np.maximum(self.log_probs, prob_1_top.T[:, :, None])
        prob_2 = prob_2.sum(axis=1).flatten()
        # prob_2 = np.sort(prob_2, axis=1)[:, remove_n:].sum(axis=1).flatten()
        # Find the mean
        # prob_2 = (log_probs + prob_1_top.T[:, :, None] * (prob_iter - 1)).sum(axis=1).flatten() / prob_iter
        prob_2_allele = np.hstack([
            # [[1],[3],[5]] -> [[1],[3],[5],[1],[3],[5]]
            np.repeat(prob_1_top_allele, self.log_probs.shape[1], axis=0),
            # [1,3,5] -> [1,3,5,1,3,5]
            np.tile(np.arange(self.log_probs.shape[1]), len(prob_1_top_allele))[:, None],
        ])
        prob_2_index = np.array(list(reversed(np.argsort(prob_2)[-self.top_n:])))

        # additional value
        prob_2_top_allele, prob_2_index = self.uniqueAllele(prob_2_allele, prob_2_index)
        prob_2_top = self.log_probs[:, prob_2_top_allele].max(axis=2)
        prob_2_value = prob_2[prob_2_index]
        prob_2_belong = np.equal(self.log_probs[:, prob_2_top_allele], prob_2_top[:, :, None])  # reads x top-n x allele_num(0/1)

        # Count of uniquely assigned
        prob_2_unique_mapped = prob_2_belong.sum(axis=2) == 1
        prob_2_unique_count = prob_2_belong.copy()
        prob_2_unique_count[np.logical_not(prob_2_unique_mapped)] = 0
        prob_2_unique_count = prob_2_unique_count.sum(axis=0)

        prob_2_fraction = (prob_2_belong / prob_2_belong.sum(axis=2)[:, :, None]).sum(axis=0)
        prob_2_fraction = prob_2_fraction / prob_2_fraction.sum(axis=1, keepdims=True)
        prob_2_fraction_norm = prob_2_fraction / self.allele_length[prob_2_top_allele]
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

        self.result.append({
            'n': len(self.result) + 1,
            'allele_id': prob_2_top_allele,
            'best_prob': prob_2_top,
            'value': prob_2_value,
            'fraction': prob_2_fraction,
            'fraction_norm': prob_2_fraction_norm,
            'unique_count': prob_2_unique_count,
            'allele_name': [list(map(self.allele_name_map_rev.get, i)) for i in prob_2_top_allele],
        })
        return self.result[-1]

    def printProb(self, n=None):
        """
        Example output:
        ```
        2
        0 -52.37587690948907
          id   1 name KIR3DS1*0130108      fraction 0.35294117647058826
          id   2 name KIR3DS1*078          fraction 0.6470588235294118
        1 -52.37587690948907
          id   0 name KIR3DS1*0130107      fraction 0.4215686274509804
          id   2 name KIR3DS1*078          fraction 0.5784313725490197
        ```
        """
        res = self.result
        if n is not None:
            res = [self.result[n]]

        for data in res:
            n = data['n']
            print(n)
            for idx, candidate in enumerate(zip(data['allele_id'],
                                                data['allele_name'],
                                                data['fraction'],
                                                data['fraction_norm'],
                                                data['unique_count'],
                                                data['value'])):
                if idx >= self.print_num:
                    break
                print(idx, candidate[-1])
                for allele_id, allele_name, frac, frac_norm, unique_count in zip(*candidate[:-1]):
                    print(f"  id {allele_id:3} name {allele_name:20s}"
                          f" fraction {frac:.5f} "
                          f" unique_count {unique_count}")

    def plot(self):
        probs = self.probs
        norm_probs = probs / probs.sum(axis=1, keepdims=True)
        # log_probs = np.log10(norm_probs)
        # plot prob
        figs = []
        figs.extend([
            px.imshow(np.log(probs), text_auto=True, aspect="auto", title=gene),
            px.imshow(norm_probs ** 0.1, text_auto=True, aspect="auto"),
            # px.imshow(np.log(norm_probs), text_auto=True, aspect="auto"),
        ])
        """
        figs.extend([dashbio.Clustergram(
            data=probs,
            column_labels=allele_names,
            hidden_labels='row',
            width=2000,
            height=1000,
            cluster="row",
        )])
        """
        return figs


class CNDist:
    def __init__(self):  # , assume_base=None):
        # const
        self.bin_num = 500
        self.max_cn = 6
        self.y0_dev = 1.5

        # parameters
        self.x_max = 1
        self.base_dev = 0.08
        # self.assume_base = assume_base
        self.assume_3DL3_diploid = True

        # result (saved for plotting)
        self.values = []
        self.loss = []
        self.base = None

    def calcProb(self, base):
        """
        Return ( CN x norm_read_count(500bin) ) array
        indicate the probility of read_count belong to the CN
        """
        x = np.linspace(0, self.x_max, self.bin_num)
        space = self.x_max / self.bin_num
        cn = np.arange(1, self.max_cn)
        y0 = norm.pdf(x, loc=0, scale=self.base_dev * self.y0_dev)
        y = np.stack([y0, *[norm.pdf(x, loc=base*n, scale=self.base_dev*n) for n in cn]])
        return y * space  # * space is to make y-sum = 1

    def fit(self, count_list):
        loss = []  # loss
        max_count = max(count_list) * 1.2

        # TODO: remove magic number
        self.base_dev *= max_count
        self.x_max = max_count 
        num_read_count, _ = np.histogram(count_list, bins=self.bin_num, range=(0, self.x_max))

        for base in np.linspace(0, self.x_max, self.bin_num):
            y = self.calcProb(base)
            max_cn_prob = y.max(axis=0)  # We will choose the highest probility of CN at that read count
            loss.append((base, np.sum( np.log(max_cn_prob + 1e-9) * num_read_count )))

        loss = np.array(loss)
        max_point = loss[np.argmax(loss[:, 1]), :]
        base = max_point[0]

        self.loss = loss
        self.base = base
        self.values = count_list

        # special case
        # if self.assume_base and base / self.assume_base > 1.7:  # TODO: maybe more CN?
        #     base /= 2
        #     self.base = base

        return base

    def assignCN(self, count_list):
        y = self.calcProb(self.base)
        max_cn_arg = y.argmax(axis=0)
        space = self.x_max / self.bin_num
        return [max_cn_arg[int(cn / space)] for cn in count_list]

    def predictCN(self, read_count):
        """ fit + assign """
        count_list = list(read_count.values())
        base = self.fit(count_list)
        cn = self.assignCN(count_list)
        return dict(zip(
            read_count.keys(),
            cn
        ))

    def predictKIRCN(self, read_count):
        """ fit + assign + 3DL3 calibrated """
        assert self.base is None
        predict_cn = self.predictCN(read_count)
        if self.assume_3DL3_diploid and predict_cn["KIR3DL3*BACKBONE"] != 2:
            print("WARNING  KIR3DL3 is not diploid, trying to fix it")
            assert predict_cn["KIR3DL3*BACKBONE"] == 1
            # we only deal with this case
            self.base /= 2

            predict_cn = dict(zip(
                read_count.keys(),
                self.assignCN(list(read_count.values()))
            ))
            assert predict_cn["KIR3DL3*BACKBONE"] == 2
        return predict_cn

    def plot(self):
        assert self.base is not None

        fig_loss = px.line(x=self.loss[:, 0], y=self.loss[:, 1])
        fig_loss.add_vline(x=self.base, line_dash="dash", line_color="gray",
                           annotation_text=f"max value={self.base}")
        fig_loss.update_layout(
            xaxis_title="read-depth",
            yaxis_title="likelihood",
        )

        fig_dist = make_subplots(specs=[[{"secondary_y": True}]])
        fig_dist.add_vline(x=self.base, line_dash="dash", line_color="gray",
                           annotation_text=f"CN=1 value={self.base}")
        fig_dist.add_trace(go.Histogram(x=list(self.values), name="Observed", nbinsx=self.bin_num, histnorm="probability"), secondary_y=True)
        # fig_dist.add_trace(go.Scatter(x=list(self.values),
        #                               y=[0 for i in range(len(self.values))],
        #                               mode='markers', name="Observed"))
        fig_dist.update_layout(
            xaxis_title="read-depth",
            yaxis_title="probility",
        )

        x = np.linspace(0, self.x_max, self.bin_num)
        y = self.calcProb(self.base)
        for n in range(len(y)):
            fig_dist.add_trace(go.Scatter(x=x, y=y[n], name=f"cn={n}"))
        return [fig_loss, fig_dist]


class KDEcut:
    def __init__(self):
        self.bandwidth = 0.05
        self.points = 100
        self.neighbor = 5  # 5 / 100 = 0.05

    def fit(self, values):
        # normalize to 0 - 1
        self.max = np.max(values)
        data = np.array(values)[:, None] / self.max
        self.kde = KernelDensity(kernel='gaussian', bandwidth=self.bandwidth).fit(data)

        # cut
        x = np.linspace(0, 1.1, self.points)
        y = self.kde.score_samples(x[:, None])
        self.local_min = x[argrelextrema(y, np.less, order=self.neighbor)[0]]
        print(self.local_min)

        # for plot
        self.data = values

    def assignCN(self, values):
        x = np.array(values) / self.max
        cn = np.searchsorted(self.local_min, x)
        return cn

    def predictCN(self, read_count):
        """ fit + assign """
        count_list = list(read_count.values())
        base = self.fit(count_list)
        cn = self.assignCN(count_list)
        return dict(zip(
            read_count.keys(),
            cn
        ))

    def plot(self):
        x = np.linspace(0, 1.1, self.points)
        y = self.kde.score_samples(x[:, None])
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig.add_trace(go.Scatter(x=x, y=y, name="kde"))
        for cn, m in enumerate(self.local_min):
            fig.add_vline(x=m, line_width=2, line_dash="dash", annotation_text=f"cn={cn}")

        fig.add_trace(go.Histogram(
            x=np.array(self.data) / self.max,
            name="Relative Depth", nbinsx=100, histnorm="probability"), secondary_y=True)
        return [fig]


class HisatKIR(Hisat2):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reference = readReference(args[0])

        # instance
        self.exon_only_cn = False
        self.cn_median = True
        self.cn_multi = False
        self.cn_cohert = False
        self.cn_dev = 0.08
        self.cn_kde = False

        # self.assume_depth = None
        self.figs = []

        self.cn_type = "sam_exon_depth"
        self.cn_type = "ans"
        self.cn_type = "variant_depth"
        self.cn_type = "sam_depth"

        self.typing_by = "likelihood"
        self.typing_by = "hisat"
        self.typing_by = "likelihood_multi"

    def plot(self):
        app = Dash(__name__)
        app.layout = html.Div([dcc.Graph(figure=f) for f in self.figs])
        app.run_server(debug=True, port=8051)

    def savePlot(self, html_filename):
        with open(html_filename, 'w') as f:
            for fig in self.figs:
                f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))

    def getGeneDepthbySamtools(self, depth_filename): # , assume_depth=None):
        df = readSamtoolsDepth(depth_filename)
        if self.exon_only_cn:
            df = cropSamtoolsDepthInExon(df, self.exon_region)

        # self.figs.append(px.line(df, x="pos", y="depth", color="gene", title="Read Depth"))
        self.figs.append(px.box(df, x="gene", y="depth", title="Read Depth"))

        return df

    def depth2CN(self, df):
        # main
        if self.cn_median:
            gene_total = df.groupby(by="gene", as_index=False)['depth'].median()
        else:
            gene_total = df.groupby(by="gene", as_index=False)['depth'].quantile(.75)
        print(gene_total)
        gene_total = dict(zip(gene_total['gene'], gene_total['depth']))

        if self.cn_kde:
            kde = KDEcut()
            gene_cn = kde.predictCN(gene_total)
            self.figs.extend(kde.plot())
        else:
            dist = CNDist()  # assume_base=self.assume_depth)
            dist.base_dev = self.cn_dev
            gene_cn = dist.predictKIRCN(gene_total)
            self.figs.extend(dist.plot())

        return gene_cn

    @staticmethod
    def selectBestTypingByRatio(result):
        # TODO: testing the new selection
        # Remove allele fraction < 0.5 / cn 
        best_id = 0
        for i in range(len(result['fraction'])):
            expect_prob = 1 / result['n']
            if all(j >= expect_prob / 2 for j in result['fraction'][i]):
                break
            best_id += 1
        if best_id >= len(result['fraction']):
            best_id = 0
        return result['allele_name'][best_id]

    def getCNFunc(self, input_gene_cn=None, return_df=False):
        if self.cn_type:
            suffix = ".cn_" + self.cn_type
        if self.cn_multi:
            suffix += "_multi"
        if not self.cn_median:
            suffix += "_p75"
        if self.cn_dev != 0.08:
            suffix += f"_dev{str(self.cn_dev).replace('.', ',')}"
        if self.cn_kde:
            suffix += f"_kde"
        if self.cn_cohert:  # this cn is pre-computed in getCNbyCohert, already saved in tsv
            suffix += "_cohert"
            return suffix, None

        if return_df:
            depth2CN = lambda i: i
        else:
            depth2CN = self.depth2CN

        # check existence ere

        if self.cn_type == "sam_depth":
            if not self.cn_multi:
                gene_cn_func = lambda name: depth2CN(self.getGeneDepthbySamtools(f"{name}.no_multi.depth.tsv"))
            else:
                gene_cn_func = lambda name: depth2CN(self.getGeneDepthbySamtools(f"{name}.depth.tsv"))

        elif self.cn_type == "sam_exon_depth":
            self.exon_only_cn = True
            if not self.cn_multi:
                gene_cn_func = lambda name: depth2CN(self.getGeneDepthbySamtools(f"{name}.no_multi.depth.tsv"))
            else:
                gene_cn_func = lambda name: depth2CN(self.getGeneDepthbySamtools(f"{name}.depth.tsv"))

        elif self.cn_type == "variant_depth":
            def tmpVariantDepth(name):
                save_reads = reconstructReadsVariants(f"{name}.id_only.json")
                if not self.cn_multi:
                    save_reads = removeMultipleAlignmentRead(save_reads)
                return depth2CN(self.getGeneDepthbySamtools(save_reads))
            gene_cn_func = tmpVariantDepth

        elif self.cn_type == "ans":
            assert input_gene_cn is not None
            gene_cn = input_gene_cn
            if "KIR2DL5*BACKBONE" in self.reference:
                gene_cn = merge2DL5CN(gene_cn)
            if "KIR2DL1S1*BACKBONE" in self.reference:
                gene_cn = merge2DL1S1CN(gene_cn)
            gene_cn_func = lambda name: {k + "*BACKBONE": v for k, v in gene_cn.items()}

        return suffix, gene_cn_func

    def plotHisatAbundance(self, hisat_report):
        fig_count = make_subplots(
            rows=1, cols=len(hisat_report),
            subplot_titles=list(map(getGeneName, hisat_report.keys())),
            shared_yaxes=True,
        )
        fig_abund = make_subplots(
            rows=1, cols=len(hisat_report),
            subplot_titles=list(map(getGeneName, hisat_report.keys())),
            shared_yaxes=True,
        )
        for i, (gene, stat) in enumerate(hisat_report.items()):
            fig_count.add_trace(
                go.Bar(y=[i[1] for i in stat['count'][:5]]),
                row=1, col=i + 1
            )
            fig_abund.add_trace(
                go.Bar(y=[i[1] for i in stat['prob'][:5]]),
                row=1, col=i + 1
            )
        
        fig_count.update_layout(title="Top-5 count",     showlegend=False)
        fig_abund.update_layout(title="Top-5 abundance", showlegend=False)
        self.figs.append(fig_count)
        self.figs.append(fig_abund)

    def getTypingFunc(self, name, gene_cn=None):
        if self.typing_by == "hisat":
            hisat_report = readReport(name + ".hisat.json")
            self.plotHisatAbundance(hisat_report)
            return ".type_hisat", lambda gene, cn: typingWithReport(hisat_report[gene], cn)

        elif self.typing_by in ["likelihood", "likelihood_multi"]:
            save_reads = reconstructReadsVariants(f"{name}.id_only.json")
            if self.typing_by == "likelihood":
                save_reads = removeMultipleAlignmentRead(save_reads)
            self.plotAlignConfusion(save_reads, no_multi=True,  title=f"Reads (gene level, no-duplication) {name}")
            self.plotAlignConfusion(save_reads,                 title=f"Reads (gene level, weighted) {name}")
            self.plotAlignConfusion(save_reads, weighted=False, title=f"Reads (gene level) {name}")
            self.plotAlignConfusion(save_reads, level="allele", title=f"Reads (allele level) {name}")

            def tmpLikelihoodTyping(gene, cn):
                reads = save_reads[gene]
                typing = AlleleTypping(reads, self.seq_len)
                result = typing.typingByNum(cn)
                alleles = self.selectBestTypingByRatio(result)
                return alleles

            return f".type_{self.typing_by}", tmpLikelihoodTyping

    def main(self, name, input_gene_cn=None, force=False):
        suffix = ".linnil1"
        new_suffix, gene_cn_func = self.getCNFunc(input_gene_cn)
        suffix += new_suffix
        print(f"{name}{suffix}.tsv")

        if os.path.exists(f"{name}{suffix}.tsv"):
            gene_cn = readCNResult(f"{name}{suffix}.tsv")
        else:
            gene_cn = gene_cn_func(name)
            df = pd.DataFrame(gene_cn.items(), columns=["gene", "cn"])
            df.to_csv(f"{name}{suffix}.tsv", sep="\t", index=False)
            self.savePlot(f"{name}{suffix}.html")
        print(gene_cn)

        new_suffix, typing_func = self.getTypingFunc(name)
        suffix += new_suffix

        if os.path.exists(f"{name}{suffix}.tsv"):
            called_alleles = readAlleleResult(f"{name}{suffix}.tsv")[0]
        else:
            called_alleles = []
            for gene, cn in gene_cn.items():
                if not cn:
                    continue
                print(gene, cn)
                # main
                called_alleles.extend(typing_func(gene, cn))
            df = pd.DataFrame([{
                'name': name,
                'alleles': "_".join(called_alleles),
            }])
            df.to_csv(f"{name}{suffix}.tsv", sep="\t", index=False)
            self.savePlot(f"{name}{suffix}.html")

        print(called_alleles)
        return suffix

    def plotAlignConfusion(self, data, level="gene", weighted=True, no_multi=False, title=""):
        gene_length = avgGeneLength(self.seq_len)

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
                                 'from_length': gene_length[gfrom] if level == "gene" else self.seq_len[gfrom],
                                 'to_length':   gene_length[gto],
                                 'value': count,
                                })
        df = pd.DataFrame(gene_from_to)
        df['norm_value_from'] = df['value'] / df['from_length']
        df['norm_value_to']   = df['value'] / df['to_length']

        order = {'from': sorted(set(df['from'])), 'to': sorted(set(df['to']))}
        self.figs.append(px.bar(df, x="from", y="value",           color="to",   text='to',   category_orders=order, title=title))
        self.figs.append(px.bar(df, x="to",   y="value",           color="from", text='from', category_orders=order,))
        self.figs.append(px.bar(df, x="from", y="norm_value_from", color="to",   text='to',   category_orders=order,))
        self.figs.append(px.bar(df, x="to",   y="norm_value_to",   color="from", text='from', category_orders=order,))
        return df

    @staticmethod
    def filterSnpOnly(alleles):
        for allele in alleles:
            if allele['typ'] != "single":
                continue
            if allele['id'].startswith("nv"):
                continue
            yield allele['id']

    def getGeneDepthbyVariant(self, data):
        # insertion deletion has more count of course
        variants = []
        for gene, reads in data.items():
            pos = []
            neg = []
            for read in reads:
                pos.extend(self.filterSnpOnly(read['lp']))
                pos.extend(self.filterSnpOnly(read['rp']))
                neg.extend(self.filterSnpOnly(read['ln']))
                neg.extend(self.filterSnpOnly(read['rn']))
            pos = Counter(pos)
            neg = Counter(neg)

            for v in (pos.keys() | neg.keys()):
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
        order_total = [i['name'] for i in sorted(variants, key=lambda i: i['total'])]
        self.figs.append( px.bar(df, x="name", y=["pos_num", "neg_num"],
                            category_orders={'name': order_total},
                            title="Depth of variants") )
        """
        self.figs.append( px.bar(df, x="name", y=["total"], color="gene",
                            category_orders={'name': order_total},
                            color_discrete_sequence=px.colors.qualitative.Dark24,
                            title="Depth of variants group by gene") )
        order_ratio = [i['name'] for i in sorted(variants, key=lambda i: i['ratio'])]
        self.figs.append( px.bar(df, x="name", y="ratio", color="gene",
                            color_discrete_sequence=px.colors.qualitative.Dark24,
                            category_orders={'name': order_ratio},
                            title="Ratio of variants") )
        """
        self.figs.append( px.box(df, x="gene", y="total",
                            category_orders={'gene': sorted(data.keys())},
                            title="Depth of each gene (no del/ins)") )

        print(df.groupby(by="gene", as_index=False).size())
        return df.rename(columns={'total': 'depth'})

    def getCNbyCohert(self, names):
        assert self.cn_cohert == False
        dfs = [self.getCNFunc(return_df=True)[1](name) for name in names]
        self.figs.append(px.histogram(pd.concat(dfs), x="depth", color="gene"))
        dfs = [df.groupby(by="gene", as_index=False)['depth'] for df in dfs] 
        if self.cn_median:
            dfs = [df.median() for df in dfs] 
        else:
            dfs = [df.quantile(.75) for df in dfs] 
        values = [v for df in dfs for v in df['depth']]
        print(values)
        self.figs.append(px.histogram(values, nbins=30))
        # print(values)

        if self.cn_kde:
            dist = KDEcut()
            dist.fit(values)
            self.figs.extend(dist.plot())
        else:
            dist = CNDist()
            dist.base_dev = self.cn_dev
            dist.fit(values)
            self.figs.extend(dist.plot())

        gene_cns = {}
        for name, df in zip(names, dfs):
            gene_cns[name] = dict(zip(
                df['gene'],
                dist.assignCN(df['depth'])
            ))
        return gene_cns

    def calcAndSaveCohertCN(self, names, output_name):
        suffix = ".linnil1"  # same in main
        new_suffix = self.getCNFunc(return_df=True)[0]
        suffix += new_suffix + "_cohert"

        if os.path.exists(f"{output_name}{suffix}.tsv"):
            gene_cn_df = pd.read_csv(f"{output_name}{suffix}.tsv", sep="\t")
        else:
            # main
            gene_cns = self.getCNbyCohert(names)
            pprint(gene_cns)

            # write all cn into tsv
            dfs = []
            for name, gene_cn in gene_cns.items():
                df = pd.DataFrame(gene_cn.items(), columns=["gene", "cn"])
                df['name'] = name
                dfs.append(df)
            gene_cn_df = pd.concat(dfs)
            gene_cn_df.to_csv(f"{output_name}{suffix}.tsv", sep="\t", index=False)
            self.savePlot(f"{output_name}{suffix}.html")

            # write cn to tsv separtely
            for name in names:
                gene_cn_df[gene_cn_df['name'] == name].to_csv(f"{name}{suffix}.tsv", sep="\t", index=False, columns=["gene", "cn"])

        return gene_cn_df


if __name__ == "__main__":
    # main("index2/kir_2100_raw.mut01", "data2/linnil1_syn_30x_seed87.00.index2_kir_2100_raw.mut01.bam")
    kir = HisatKIR("index/kir_2100_raw.mut01")

    cohert_name = "data/linnil1_syn_30x_seed87"
    name = cohert_name.split("/")[1]
    ans = EvaluateKIR(f"{name}/{name}.summary.csv")

    n = 10
    names = [f"data/linnil1_syn_30x_seed87.{i:02d}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr" for i in range(n)]
    gene_cns = kir.calcAndSaveCohertCN(names, "tmp")
    kir.plot()
    exit()

    for i in range(n):
        id = f"{i:02d}"
        name = f"data/linnil1_syn_30x_seed87.{id}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr"
        gene_cn = dict(Counter(map(getGeneName, ans.getAns(id))))
        suffix = kir.main(name, input_gene_cn=gene_cn)

    called_alleles_dict = {}
    predict_list = []
    for i in range(n):
        id = f"{i:02d}"
        name = f"data/linnil1_syn_30x_seed87.{id}.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr"
        predict = readAlleleResult(name + suffix + ".tsv")[0]
        predict['id'] = id
        predict_list.append(predict)
        called_alleles_dict[id] = predict['alleles']

    df = pd.DataFrame(predict_list)
    # df.to_csv(f"{cohert_name}_merge{suffix}.{kir.typing_by}.tsv", index=False, sep="\t")

    ans.compareCohert(called_alleles_dict, skip_empty=True)
    kir.plot()
