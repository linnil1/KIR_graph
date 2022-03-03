import pandas as pd


def evaluteList(a, b):
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


def evalute(answer_list, predict_list):
    ans_list = sorted(answer_list)
    predict_list = sorted(predict_list)
    match_pair = []

    for gene in set([i.split('*')[0] for i in ans_list + predict_list]):
        match_pair.extend(evaluteList(
            [i for i in ans_list if i.startswith(gene)],
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


def getSummaryAnswer(id):
    data = pd.read_csv("linnil1_syn_full/summary.csv", sep="\t", header=None)
    data = list(data[2].str.split("_"))[0]
    print(data)
    return data


def getPingResult(id):
    data = pd.read_csv("/home/linnil1/kir/PING/PING_test/linnil1_syn_full_result/finalAlleleCalls.csv", sep=',')
    print(data)
    pred = []
    for i in data.loc[id - 1]:
        if "linnil1" not in i:
            pred.extend(i.split('+'))
    return pred


id = 2
evalute(getSummaryAnswer(id), getPingResult(id))
