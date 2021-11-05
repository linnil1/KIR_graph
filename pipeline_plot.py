import asyncio
import subprocess
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor
import re
from matplotlib.ticker import MaxNLocator
import numpy as np
from collections import defaultdict
from Bio import SeqIO


def getPairNum(name):
    cmd = f"docker run -it --rm --name {name.replace('/', '_')} -v $PWD:/app -w /app samtools samtools flagstat {name}"
    print(name, cmd)
    a = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    # print(a.stdout)
    num_pair = 0
    num_total = 0
    for i in a.stdout.decode().split("\n"):
        # 122050 + 0 properly paired (100.00% : N/A)
        # 122050 + 0 in total (QC-passed reads + QC-failed reads)
        # 0 + 0 secondary
        if "in total" in i:
            num_total = int(i.split()[0])
        if "secondary" in i:
            num_second = int(i.split()[0])
        if "properly paired" in i:
            num_pair = int(i.split()[0])
    return num_total - num_second, num_second, num_pair


def plot_bam_mapping():
    perc_secd = []
    perc_pair = []
    n = 10

    exc = []
    with ThreadPoolExecutor(max_workers=50) as executor:
        # tot, sec, pair = getPairNum(name)
        for i in range(n):
            name = f"data/linnil1_syn.{i:02d}.merge.sam"
            exc.append(executor.submit(getPairNum, name))
            # exc.append(getPairNum(name))

            name = f"data/linnil1_syn.{i:02d}.split.sam"
            exc.append(executor.submit(getPairNum, name)) # exc.append(getPairNum(name))

            name = f"data/linnil1_syn.{i:02d}.linear.sam"
            exc.append(executor.submit(getPairNum, name))
            # exc.append(getPairNum(name))

            name = f"data/linnil1_syn.{i:02d}.full.sam"
            exc.append(executor.submit(getPairNum, name))
            # exc.append(getPairNum(name))

        for job in exc:
            # tot, sec, pair = getPairNum(name)
            tot, sec, pair = job.result()
            # tot, sec, pair = job
            perc_pair.append(pair / tot)
            perc_secd.append(sec / tot)


    # plt.subplot(1, 3, 1)
    plt.title("Proper Paired percetage")
    plt.boxplot([perc_pair[0::4], perc_pair[1::4], perc_pair[2::4], perc_pair[3::4]])
    plt.gca().set_xticklabels(['merge', 'split', 'linear', 'full'])
    plt.show()

    # plt.subplot(1, 3, 2)
    plt.title("Secondary percetage")
    plt.boxplot([perc_secd[0::4], perc_secd[1::4], perc_secd[2::4], perc_secd[3::4]])
    plt.gca().set_xticklabels(['merge', 'split', 'linear', 'full'])
    plt.show()


def plot_full_acc():
    arr_acc = []
    arr_acc_gene = []
    for i in range(10):
        tot, acc, acc_gene = 0, 0, 0
        for i in open(f"data/linnil1_syn.{i:02d}.full.sam"):
            if i and "@" == i[0]:
                continue
            row = i.split()
            if len(row) < 2:
                print(i)
            ans = row[0].split('-')[0]
            pred = row[2]
            tot += 1
            if ans == pred:
                acc += 1
            if ans.split("*")[0] == pred.split("*")[0]:
                acc_gene += 1
        arr_acc.append(acc / tot)
        arr_acc_gene.append(acc_gene / tot)

    # plt.subplot(1, 3, 3)
    plt.title("Accuracy for mapped on all kir star alleles")
    plt.boxplot([arr_acc, arr_acc_gene])
    plt.gca().set_xticklabels(['Allele', 'Gene'])
    plt.show()


def plot_each_iter():
    start_save_count = False
    loss = []
    tpm  = []
    name = []
    rec = []
    acc = []
    for i in open("tmp_out1"):
        if "Rank" in i:
            if rec:
                if len(rec) > 28:
                    # rec, tpm, name = list(zip(*sorted(zip(rec, tpm, name))[::-1]))
                    plt.suptitle(f"allele num = {len(rec)}")
                    plt.subplot(1,2,1)
                    plt.title(f"perpotion")
                    plt.bar(range(len(rec)), rec)
                    ax = plt.gca()
                    ax.set_xticks(range(len(rec)))
                    ax.set_xticklabels(name, rotation=90)
                    plt.subplot(1,2,2)
                    plt.title(f"tpm")
                    plt.bar(range(len(rec)), tpm)
                    ax = plt.gca()
                    ax.set_xticks(range(len(rec)))
                    ax.set_xticklabels(name, rotation=90)
                    plt.show()
            start_save_count = False
            rec = []
            tpm = []
            name = []
        if i.strip() and start_save_count and "Allele" in i:
            name.append(i.split()[1])
            rec.append(float(i.split()[-1]))
            tpm.append(float(i.split()[-3]))
        if i.strip() and start_save_count and "ACC" in i:
            acc.append(float(i.split()[1]))
        if "Rank 0" in i:
            loss.append(float(i.split()[-1]))
            start_save_count = True


    plt.grid()
    x = np.arange(len(loss)) + 1
    plt.xticks(x)
    plt.plot(x, loss, 'b.-')
    plt.ylabel("loss", color='b')
    ax2 = plt.gca().twinx()
    ax2.plot(x, acc, 'g.-')
    ax2.set_ylabel("Acc", color='g')
    plt.show()



def read_full_pair(fname, strict=False):
    read_map = defaultdict(list)
    for line in open(fname):
        if line and "@" == line[0]:
            continue
        if strict and "150M" not in line:
            continue
        row = line.split()
        ans = row[0].split('-')[0]
        if int(row[1]) & 64:
            read_map[row[0] + "/1"].append(row[2])
        else:
            read_map[row[0] + "/2"].append(row[2])

    pair_map = {}
    for i in set(read_map.keys()):
        if i[-1] == "1":
            id = i[:-2]
            pair_map[id] = set(read_map[id + "/1"]) & set(read_map[id + "/2"])
    return pair_map


def plotSplitAcc():
    i = 0
    pair_map = read_full_pair(f"data/linnil1_syn.{i:02d}.full.sam")
    pair_num = [len(i) for i in pair_map.values()]
    plt.title("Number of multiple alignment")
    plt.boxplot(pair_num)
    plt.show()


def plotAnsDepth():
    seqs = SeqIO.parse("kir_merge_sequences.fa", "fasta")
    seq_len = {seq.id: len(seq.seq) for seq in seqs}

    for i in range(10):
        allele_count = defaultdict(int)
        for line in open(f"linnil1_syn/linnil1_syn.{i:02d}.read..sam"):
            if line[0] == "@":
                continue
            row = line.split()
            allele_count[row[0].split("-")[0]] += 1

        alleles = [(a, b, b / seq_len[a] * 150) for a, b in allele_count.items()]
        plotDepth(sorted(alleles, key=lambda j: -j[2]))


def plotAnsDepthWithMulti():
    seqs = SeqIO.parse("kir_merge_sequences.fa", "fasta")
    seq_len = {seq.id: len(seq.seq) for seq in seqs}

    alleles_set = []
    i = 0
    for line in open(f"linnil1_syn/linnil1_syn.{i:02d}.read..sam"):
        if line[:3] == "@SQ":
            alleles_set.append(line.split("\t")[1][3:])
    alleles_set = set(alleles_set)

    allele_count = defaultdict(int)
    pair_map = read_full_pair(f"data/linnil1_syn.{i:02d}.full.sam", strict=True)
    for id, alleles in pair_map.items():
        ans_alleles = alleles_set & alleles
        for i in ans_alleles:
            allele_count[i] += 1 / len(ans_alleles)

    alleles = [(a, b, b / seq_len[a] * 150) for a, b in allele_count.items()]
    plotDepth(sorted(alleles, key=lambda j: -j[2]))


def plotDepth(alleles, need_sort=False):
    # alleles = [(name, perpotion, tpm),]
    if need_sort:
        alleles = sorted(alleles, key=lambda i: -i[2])
    n = list(range(len(alleles)))
    plt.suptitle(f"Allele num = {len(n)}")
    allele_names, allele_count, allele_tpm = zip(*alleles)
    plt.subplot(1,2,1)
    plt.title(f"Depth")
    plt.bar(n, allele_tpm)
    ax = plt.gca()
    ax.set_xticks(n)
    ax.set_xticklabels(allele_names, rotation=90)
    plt.subplot(1,2,2)
    plt.title(f"perpotion")
    plt.bar(n, allele_count)
    ax = plt.gca()
    ax.set_xticks(n)
    ax.set_xticklabels(allele_names, rotation=90)
    plt.show()


if __name__ == "__main__":
    # plot_bam_mapping()
    # plot_full_acc()
    plotAnsDepth() 
    # plotAnsDepthWithMulti() 
