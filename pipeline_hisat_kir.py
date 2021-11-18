import numpy as np
import json
from collections import defaultdict
from pprint import pprint
import matplotlib.pyplot as plt
import glob
import sys
import os
import re
from Bio import SeqIO
import pandas as pd
from pipeline_plot import plotDepth


# Per index per sample per gene
class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)


class HisatTyping:
    def __init__(self, index="kir_merge"):
        # sequences length (for normalize)
        self.index = index
        for seq in SeqIO.parse(f"{self.index}_sequences.fa", "fasta")
            self.seq_len[seq.id] = len(seq.seq)
        for seq in SeqIO.parse(f"{self.index}_backbone.fa", "fasta"):
            self.seq_len[seq.id] = len(seq.seq)

    def calcProbPerPair(self):
        print(f"Read {self.name}.{self.gene}.json")
        # load data
        allele_count_per_read = json.load(open(f'{self.name}.{self.gene}.json'))

        """
        # separate mode
        if not len(allele_count_per_read):
            return
        if len(allele_count_per_read) < 100:  # min reads
            return
        """
        # name -> index
        id_map = set()
        for j in allele_count_per_read:
            id_map = id_map | dict(j['p']).keys() | dict(j['n']).keys()
        id_map = {i: ind for ind, i in enumerate(id_map)}

        # 0.999 when correct, 0.001 when incorrect
        allele_prob_per_read = []
        for j in allele_count_per_read:
            p, n = dict(j['p']), dict(j['n'])
            arr = np.ones(len(id_map))
            for i in p:
                arr[id_map[i]] *= 0.999 ** p[i]
            for i in n:
                arr[id_map[i]] *= 0.001 ** n[i]
            allele_prob_per_read.append(arr)
        print("Per seq probility Done")

        # save it
        self.id_map = id_map
        self.allele_prob_per_read = np.array(allele_prob_per_read)
        np.savez(f'{self.name}.{self.gene}.tmpprob.npz',  read=self.allele_prob_per_read, idmap=self.id_map)

    def calcProbAlleles(self):
        print(f"Calculate the probility of all reads for all alleles")
        allele_prob_per_read = self.allele_prob_per_read
        id_map = self.id_map
        top_n = 10    # reduce computation
        save_n = 100  # reduce data
        if "merge" in self.index:
            max_iter = 40
        elif "split" in self.index:
            max_iter = 5
        self.allele_prob_n_iter = []

        # one allele
        allele_prob_all_1 = np.log(allele_prob_per_read + 1e-100).sum(axis=0) / len(allele_prob_per_read)
        allele_prob_all_1_list = [(allele_prob_all_1[v1], [k1]) for k1, v1 in id_map.items() if not k1.endswith("BACKBONE")]
        allele_prob_all_1_list = list(map(list, sorted(allele_prob_all_1_list, key=lambda i:-i[0])))
        # pprint(allele_prob_all_1_list[:3])
        for max_allele in allele_prob_all_1_list[:top_n]:
            max_allele.append([np.mean(allele_prob_per_read.max(axis=1) == \
                                       allele_prob_per_read[:, id_map[max_allele[1][0]]]),
                               np.sum (allele_prob_per_read.max(axis=1) == \
                                       allele_prob_per_read[:, id_map[max_allele[1][0]]])])
            print(max_allele)
        self.allele_prob_n_iter.append(allele_prob_all_1_list[:save_n])

        # more than one allele
        """
        self.allele_prob_n_iter = [
            [
            # my
            [0, ['KIR2DL1*001', 'KIR2DL1*007', 'KIR2DL1*041', 'KIR2DL2*003', 'KIR2DL3*006', 'KIR2DL4*008', 'KIR2DL4*011', 'KIR2DL4*049', 'KIR2DL4*052', 'KIR2DL5A*031', 'KIR2DL5B*004', 'KIR2DL5B*007', 'KIR2DL5B*008', 'KIR2DP1*002', 'KIR2DP1*007', 'KIR2DP1*009', 'KIR2DS1*012', 'KIR2DS2*019', 'KIR2DS3*010', 'KIR2DS3*015', 'KIR2DS4*004', 'KIR2DS4*006', 'KIR2DS5*037', 'KIR3DL2*010', 'KIR3DL2*018', 'KIR3DL3*009', 'KIR3DL3*025', 'KIR3DL3*048', 'KIR3DP1*015', 'KIR3DS1*013']],
            # real
            [1, ['KIR2DL1*001', 'KIR2DL1*007', 'KIR2DL1*055', 'KIR2DL2*003', 'KIR2DL3*006', 'KIR2DL4*005', 'KIR2DL4*008', 'KIR2DL4*011', 'KIR2DL5A*012', 'KIR2DL5A*029', 'KIR2DL5A*033', 'KIR2DL5B*004', 'KIR2DP1*002', 'KIR2DP1*008', 'KIR2DP1*009', 'KIR2DS1*011', 'KIR2DS1*012', 'KIR2DS2*019', 'KIR2DS3*010', 'KIR2DS3*015', 'KIR2DS5*020', 'KIR2DS5*022', 'KIR3DL2*010', 'KIR3DL2*018', 'KIR3DL3*001', 'KIR3DL3*041', 'KIR3DP1*022', 'KIR3DP1*041', 'KIR3DP1*044', 'KIR3DS1*013']], # , 'KIR3DS1*078']
            ]
        ]
        for iter_i in range(1):
        """
        for iter_i in range(max_iter):
            print(f"Iter {iter_i}")
            allele_prob_all_prev_list = self.allele_prob_n_iter[-1]
            first_n = min(10, len(allele_prob_all_prev_list))
            # max
            allele_prob_all_now = np.array([
                np.max([allele_prob_per_read[:, id_map[j]] for j in i[1]], axis=0)
                            for i in allele_prob_all_prev_list[:first_n]]).T
            allele_prob_all_now = np.max([
                np.tile(allele_prob_per_read[:, None, :], [first_n,   1]),
                np.tile(allele_prob_all_now[:, :, None   ], [1, len(id_map)])], axis=0)
            """
            # sum
            allele_prob_all_now = np.array([
                np.sum([allele_prob_per_read[:, id_map[j]] for j in i[1]], axis=0)
                            for i in allele_prob_all_prev_list[:first_n]]).T
            allele_prob_all_now = np.sum([
                np.tile(allele_prob_per_read[:, None, :], [first_n,   1]),
                np.tile(allele_prob_all_now[:, :, None   ], [1, len(id_map)])], axis=0)
            """
            # all
            allele_prob_all_now = np.log(allele_prob_all_now + 1e-100).sum(axis=0) \
                                / len(allele_prob_per_read)
            # Note: ignore duplicated
            allele_prob_all_now_list = [
                    (allele_prob_all_now[v1, v2],
                     sorted([*allele_prob_all_prev_list[v1][1], k3])) \
                        for v1 in range(first_n) \
                            for k3, v2 in id_map.items() \
                                if not k3.endswith("BACKBONE") and k3 not in allele_prob_all_prev_list[v1][1]]
            allele_prob_all_now_list = list(map(list,
                sorted(allele_prob_all_now_list, key=lambda i:-i[0])))
            pprint(allele_prob_all_now_list[:2])
            self.allele_prob_n_iter.append(allele_prob_all_now_list[:save_n])

        # save
        json.dump(self.allele_prob_n_iter, open(f'{self.name}.{self.gene}.em.json', 'w'), cls=MyEncoder)

    def allelesAssign(self, alleles, normalize=True, map_name=False):
        """
        Input:
          ["KIR2DL1*001", "KIR2DL1*002"]

        Output:
          Basic: [ [0, 1], [1,0], [1, 1] ]
          Normalize: [ [0, 1], [1,0], [0.5, 0.5] ]
          MapName: [ ["002"], ["001"], ["002", "001"] ]
        """
        ids = [self.id_map[i] for i in alleles]
        reads_prob = self.allele_prob_per_read[:, ids]
        max_pos = np.equal(reads_prob, np.max(reads_prob, axis=1)[:, None])
        if map_name:
            ids = np.array(alleles)
            names = [ids[i == 1] for i in max_pos]
            return names

        if normalize:
            max_pos = max_pos / max_pos.sum(axis=1)[:, None]
        return max_pos

    def cutProbThreshold(self):
        print(f"Cut the threshold")
        allele_prob_per_read = self.allele_prob_per_read
        prev_loss = -np.inf
        prev_tpm = []
        prev_allele_list = []
        for allele_prob_all_now_list in self.allele_prob_n_iter:
            # rank 0
            tmp = allele_prob_all_now_list[0]
            loss, alleles = tmp[0], tmp[1]
            # count tpm
            count = self.allelesAssign(alleles).sum(axis=0)
            read_len = 150 * 2
            alleles_abundance = list(zip(
                # allele name
                alleles,
                # perpotion
                map(lambda i: count[i] / sum(count), range(len(alleles))),
                # depth
                map(lambda i: count[i] * read_len / self.seq_len[alleles[i]], range(len(alleles))),
            ))
            tpm = [i[2] for i in alleles_abundance]

            # main IF
            # log(likelihood) -> All < 0
            # larger better e.g. -3 > -5
            if prev_loss > loss * 1.01 and np.min(tpm / np.median(tpm)) < 0.6:
                print(f"{len(alleles)} Fail")
                pprint(alleles_abundance)
                break
            else:
                print(f"{len(alleles)} OK (medium={np.median(tpm)})")
                # ok
                prev_loss = loss
                prev_tpm = tpm
                prev_allele_list = alleles

        pred_allele = self.allelesAssign(prev_allele_list, map_name=True)
        print("size=", len(prev_allele_list))
        print("alleles=", prev_allele_list)
        summary = []
        for i in range(len(prev_tpm)):
            this_allele_cn = round(prev_tpm[i] / np.median(prev_tpm))
            print("allele =", prev_allele_list[i],
                  "tpm =", prev_tpm[i],
                  "tpm(norm) =", prev_tpm[i] / np.median(prev_tpm),
                  "CN =", this_allele_cn)
            for _ in range(this_allele_cn):
                summary.append(prev_allele_list[i])
        print("tot=", len(summary))
        return summary, [np.random.choice(i) for i in pred_allele]

    def printRank(self):
        allele_prob_per_read = self.allele_prob_per_read
        for allele_prob_all_now_list in self.allele_prob_n_iter:
            # print info
            for i in range(min(2, len(allele_prob_all_now_list))):
                allele = allele_prob_all_now_list[i]
                print(f"Size: {len(allele[1])} Rank {i} Loss {allele[0]}")

                # reduce
                max_pos = self.allelesAssign(allele[1])
                count = max_pos.sum(axis=0)
                for tpm, num, allele_name in sorted(zip(
                        map(lambda i: count[i] / self.seq_len[allele[1][i]], range(len(allele[1]))),
                        map(lambda i: count[i] / sum(count), range(len(allele[1]))),
                        allele[1]))[::-1]:
                    print(f"  Allele: {allele_name} TPM: {tpm} Propotion {num}")

    def likehoodPlot(self):
        num = []
        loss = []
        for allele_iter in self.allele_prob_n_iter:
            allele_list = allele_iter[0]
            loss.append(allele_list[0])
            num.append(len(allele_list[1]))
        print(np.array(loss))
        plt.plot(num, loss, '.-', label="linnil1")
        plt.ylabel("log likelihood")
        plt.xlabel("Number of allele")
        plt.legend()
        plt.show()

    def likehoodReport(self):
        allele_list = []
        for i in open("data/linnil1_syn_full.00.merge.pair.report"):
            if "ranked" in i:
                allele_list.append(i.split()[2])

        loss = []
        for i in range(1, len(allele_list)):
            loss.append(
                np.log(1e-100 + \
                       self.allele_prob_per_read[:, [self.id_map[j] for j in allele_list[:i]]]
                           .sum(axis=1))
                  .mean())
        print(np.array(loss))
        plt.plot(range(1, len(allele_list)), loss, '.-', label="report")

    def alleleDistPlot(self):
        arr_dists = []

        arr_multialigns = []
        for allele_prob_all_now_list in self.allele_prob_n_iter:
            allele = allele_prob_all_now_list[0]
            max_pos = self.allelesAssign(allele[1])
            count = max_pos.sum(axis=0)

            read_len = 150 * 2
            alleles = list(zip(
                # allele name
                allele[1],
                # perpotion
                map(lambda i: count[i] / sum(count), range(len(allele[1]))),
                # depth
                map(lambda i: count[i] * read_len / self.seq_len[allele[1][i]], range(len(allele[1]))),
            ))
            arr_dists.append(alleles)
            arr_multialigns.append([len(i) for i in self.allelesAssign(allele[1], map_name=True)])

        plt.title("Number of Multple Allele Assignment per reads")
        plt.boxplot(arr_multialigns, showfliers=False)
        plt.show()

        # plot distructibon per iteration
        json.dump(arr_dists, open(f'tmp.json', 'w'))
        num = 1
        for alleles in arr_dists[30-8:30+8]:
            alleles = sorted(alleles, key=lambda i: -i[2])
            plt.subplot(4, 4, num)
            plt.title(len(alleles))
            plt.bar(range(len(alleles)), [i[2] for i in alleles])
            num += 1
        plt.show()
        return

        # plot specific distribution
        for dist in arr_dists[29:33]:
            print(dist)
            print("len = ", len(dist))
            for name, prop, dep in dist:
                print(f"  Allele: {name} Depth: {dep} Propotion {prop}")
            plotDepth(dist, need_sort=True)

    def mainPerGene(self):
        print(f"Typing for {self.gene} in {self.name}")
        # self.evaluateHisatMap()

        if os.path.exists(f'{self.name}.{self.gene}.tmpprob.npz'):
            v = np.load(f'{self.name}.{self.gene}.tmpprob.npz', allow_pickle=True)
            self.allele_prob_per_read = np.array(v['read'])
            self.id_map = json.loads(v['idmap'].__str__().replace("'", '"'))
        else:
            self.calcProbPerPair()

        if os.path.exists(f'{self.name}.{self.gene}.em.json'):
            self.allele_prob_n_iter = json.load(open(f'{self.name}.{self.gene}.em.json'))
        else:
            self.calcProbAlleles()
        # self.evaluateByName()
        # self.likehoodReport()  # run this beofre likehoodPlot
        # self.likehoodPlot()
        # self.alleleDistPlot()
        called_allele, assign_allele = self.cutProbThreshold()

        # save to summary allele
        self.summary_allele.extend(called_allele)

        # save to bam tmp array
        for i in set(called_allele):
            self.bam_headers.append(f"@RG\tID:{i}\n")
        if len(assign_allele) > 0:
            assert len(self.bam_records_per_gene[self.gene]) == 2 * len(assign_allele)
        for i, allele_name in enumerate(assign_allele):
            self.bam_new.append(self.bam_records_per_gene[self.gene][i * 2    ] + "\tRG:Z:" + allele_name + "\n")
            self.bam_new.append(self.bam_records_per_gene[self.gene][i * 2 + 1] + "\tRG:Z:" + allele_name + "\n")

    def readCounts(self):
        num_read_per_gene = {}
        for gene, reads in self.bam_records_per_gene.items():
            num_read_per_gene[gene] = len(reads)
        num_read_per_gene = sorted(num_read_per_gene.items(), key=lambda i: -i[1])
        pprint(num_read_per_gene)

        plt.subplot(1, 2, 1)
        plt.title("Read count")
        plt.bar(range(len(num_read_per_gene)), [i[1] for i in num_read_per_gene])
        ax = plt.gca()
        ax.set_xticks(range(len(num_read_per_gene)))
        ax.set_xticklabels([i[0] for i in num_read_per_gene], rotation=90)

        num_read_per_gene = {}
        for gene, reads in self.bam_records_per_gene.items():
            num_read_per_gene[gene] = len(reads) / self.seq_len[gene + "*BACKBONE"]
        num_read_per_gene = sorted(num_read_per_gene.items(), key=lambda i: -i[1])
        pprint(num_read_per_gene)

        plt.subplot(1, 2, 2)
        plt.title("Read count / length")
        plt.bar(range(len(num_read_per_gene)), [i[1] for i in num_read_per_gene])
        ax = plt.gca()
        ax.set_xticks(range(len(num_read_per_gene)))
        ax.set_xticklabels([i[0] for i in num_read_per_gene], rotation=90)
        plt.show()


    def mainPerSample(self, name):
        self.name = name

        # split mode
        if "split" in self.index:
            genes  = ["KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B",
                      "KIR2DP1", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4",
                      "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DP1",
                      "KIR3DS1"]
        # merge mode
        elif "merge" in self.index:
            genes  = ["KIR"]
        else:
            raise "Not in split or merge mode"

        # reads bam
        self.readBam()
        self.summary_allele = []

        self.readCounts()  # split mode
        exit()

        # main
        for gene in genes:
            self.gene = gene
            self.mainPerGene()

        # tmp
        self.writeBam()
        self.writeSummary()

    def readBam(self):
        self.bam_records_per_gene = defaultdict(list)
        self.bam_headers = []
        self.bam_new = []
        f_sam = open(f"{self.name}.sam")
        for i in f_sam:
            if i[0] == "@":
                self.bam_headers.append(i)
            else:
                gene = re.findall("\t(KIR.*?)\*BACKBONE\t", i)[0]
                self.bam_records_per_gene[gene].append(i.strip())

    def writeBam(self):
        # Write all things into one bam
        newname = f"{self.name}.group"
        print(f"Write to {newname}.bam")
        f_sam = open(f"{newname}.sam", "w")
        f_sam.writelines(self.bam_headers)
        f_sam.writelines(self.bam_new)
        f_sam.close()
        os.system(f"samtools sort -@30 {newname}.sam -o {newname}.bam")
        os.system(f"samtools index {newname}.bam")

    def writeSummary(self):
        # Write all things into one bam
        newname = f"{self.name}.group.txt"
        print(f"Write to {newname}")
        with open(newname, "w") as f:
            f.write(f"{self.name}," + ",".join(self.summary_allele) + "\n")

    def evaluateByName(self):
        # data = pd.read_csv("ping_syn_data/synSeq.info.txt", sep="\t", header=None)
        for allele_cands in self.allele_prob_n_iter[20:35]:
            allele = allele_cands[0]
            pos_names = self.allelesAssign(allele[1], map_name=True)
            tot, acc, acc_gene = len(pos_names), 0, 0
            for i in range(len(pos_names)):
                if self.bam_records_per_gene[self.gene][i * 2].split("-")[0] in pos_names[i]:
                    acc += 1
                if self.bam_records_per_gene[self.gene][i * 2].split("*")[0] in set(
                        [j.split("*")[0] for j in pos_names[i]]):
                    acc_gene += 1
            print("Alleles", len(allele[1]))
            print("Allele ACC:", acc / tot)
            print("Gene ACC:", acc_gene / tot)

    def evalute(self, id, predict_list):
        data = pd.read_csv("linnil1_syn_full/summary.csv", sep="\t", header=None)
        data = list(data[2].str.split("_"))
        ans_list = sorted(data[id])
        predict_list = sorted(predict_list)
        match_pair = []

        for gene in set([i.split('*')[0] for i in ans_list + predict_list]):
            match_pair.extend(
                    self.evaluteList([i for i in ans_list if i.startswith(gene)],
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

    def evaluteList(self, a, b):
        match_pair = []
        i = 0
        for allele in list(b):
            if allele in a:
                a.remove(allele)
                b.remove(allele)
                match_pair.append(("Match", allele, allele))

        for allele_b in list(b):
            for allele_a in a:
                if allele_a[:allele_a.index("*") + 3] == allele_b[:allele_b.index("*") + 3]:
                    a.remove(allele_a)
                    b.remove(allele_b)
                    match_pair.append(("Mismatch", allele_a, allele_b))
                    break

        for allele in a:
            match_pair.append(("FN", allele, None))
        for allele in b:
            match_pair.append(("FP", None, allele))
        return match_pair

    def evaluateHisatMap(self):
        print("Read mapping...")
        allele_count_per_read = json.load(open(f'{self.name}.{self.gene}.json'))
        tot = 0
        posneg_and_acc = 0
        posneg_max_acc = 0
        posneg_maxsec_acc = 0

        print_seq = True
        if print_seq:
            from pyHLAMSA import Genemsa
            msa = Genemsa.load_msa(f"{self.index}.save.fa", f"{self.index}.save.gff")
            ans_sam = open("linnil1_syn/linnil1_syn_full.00.read..sam")
            # ans_pos = defaultdict(list)
                    # ans_pos[rec[0]].append(int(rec[3]))
            ans_pos = {}
            ans_pos_rev = {}
            for line in ans_sam:
                if line.startswith("@"):
                    continue
                rec = line.split()
                if (int(rec[1]) & 16) != 0:
                    ans_pos[rec[0]] = int(rec[3])
                else:
                    ans_pos_rev[rec[0]] = int(rec[3])

        for i in range(len(allele_count_per_read)):
            tot += 1
            p, n = allele_count_per_read[i]['p'], allele_count_per_read[i]['n']
            allele_and = set([i[0] for i in p]) - set([i[0] for i in n])
            allele_sum = dict(p)
            for name, count in n:
                if name not in allele_sum:
                    allele_sum[name] = -count
                else:
                    allele_sum[name] += -count
            if len(p):
                max_count = max(allele_sum.values())
                allele_sum_sec = set([name for name, count in allele_sum.items() if count >= max_count - 1])
                allele_sum = set([name for name, count in allele_sum.items() if count == max_count])
            else:
                allele_sum = set()
                allele_sum_sec = set()
            ans_allele = self.bam_records_per_gene[self.gene][i * 2].split('-')[0]
            if ans_allele in allele_and:
                posneg_and_acc += 1
            if ans_allele in allele_sum:
                posneg_max_acc += 1
            else:
                # save in tmp4
                print(self.bam_records_per_gene[self.gene][i * 2])
                print(self.bam_records_per_gene[self.gene][i * 2 + 1])
                print("Predict", allele_sum)

                if print_seq:
                    # real length to msa length
                    count = 0
                    pos_map = {}
                    pos_map_back = {}
                    for msa_pos, c in enumerate(msa.get(ans_allele)):
                        if c != "-":
                            pos_map[count] = msa_pos
                            pos_map_back[msa_pos] = count
                            count += 1
                        else:
                            pos_map_back[msa_pos] = count

                    allele_sum.add(ans_allele)
                    select_msa = msa.select_allele("|".join(allele_sum).replace("*", r"\*"))

                    rec = self.bam_records_per_gene[self.gene][i * 2].split()
                    rec_rev = self.bam_records_per_gene[self.gene][i * 2 + 1].split()
                    if (int(rec[1]) & 16) == 0:
                        rec, rec_rev = rec_rev, rec

                    # + strand
                    print(" Read+              ", end="")
                    for j in range(0, len(rec[9]), 10):
                        print(rec[9][j:j+10], end=' ')
                    print("")

                    pos = int(rec[3]) - 1
                    pos1 = pos_map[min(pos_map_back[pos] + 149, count - 1)] + 1
                    print(" Pos                ", pos + 1, pos1)  # included
                    print(select_msa[pos:pos1].format_alignment_diff(ans_allele))
                    if pos != pos_map[ans_pos[rec[0]]] - 1:
                        pos = pos_map[ans_pos[rec[0]]] - 1
                        pos1 = pos_map[min(ans_pos[rec[0]] + 148, count - 1)] + 1
                        print(" Real Map           ")
                        print(" Pos                ", pos + 1, pos1)  # included
                        print(select_msa[pos:pos1].format_alignment_diff(ans_allele))

                    # - strand
                    print(" Read+              ", end="")
                    for j in range(0, len(rec_rev[9]), 10):
                        print(rec_rev[9][j:j+10], end=' ')
                    print("")

                    pos = int(rec_rev[3]) - 1
                    pos1 = pos_map[min(pos_map_back[pos] + 149, count - 1)] + 1
                    print(" Pos                ", pos + 1, pos1)  # included
                    print(select_msa[pos:pos1].format_alignment_diff(ans_allele))

                    break

            if ans_allele in allele_sum_sec:
                posneg_maxsec_acc += 1

        posneg_and_acc /= tot
        posneg_max_acc /= tot
        print("Read Mapping Allele Acc(AND)", posneg_and_acc)
        print("Read Mapping Allele Acc(SUM)", posneg_max_acc)
        print("Read Mapping Allele Acc(SUM)(allow 1 error)", posneg_maxsec_acc / tot)
        return posneg_and_acc, posneg_max_acc


if __name__ == "__main__":
    # names = list(set(map(lambda i: i.split(".KIR")[0], glob.glob("data/synSeq.hiseq.dp50.rl150.*pair.tmp.*.json"))))
    # name = f"data/synSeq.hiseq.dp50.rl150.1.split.pair.tmp"
    # name = sys.argv[1]
    # name = f"data/synSeq.hiseq.dp50.rl150.1.merge.pair.tmp"

    i = 0
    names = []
    for i in range(10):
        # names.append(f"data/linnil1_syn.0{i}.merge.pair.tmp")
        # names.append(f"data/linnil1_syn_full.0{i}.merge.pair.tmp")
        names.append(f"data/linnil1_syn_full.0{i}.split.pair.tmp")

    names = names[:1]
    typ = HisatTyping("kir_merge_full")
    # typ = HisatTyping("kir_split_full")
    for name in names:
        typ.mainPerSample(name)

    '''
    with open("tmp.summary.txt", "w") as f:
        summarys = [open(f"{name}.group.txt") for name in names]
        f.writelines(summarys)
    '''
