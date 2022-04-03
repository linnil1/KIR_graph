import pandas as pd
from pprint import pprint


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


def evalute(answer_list, predict_list, field=5):
    ans_list = sorted(answer_list)
    predict_list = sorted(predict_list)
    match_pair = []

    # evaluate per same KIR gene
    for gene in set(i.split('*')[0] for i in ans_list + predict_list):
        match_pair.extend(evaluteList(
            [i[:min(len(i), i.index("*")+1+field)] for i in ans_list if i.startswith(gene)],
            [i[:min(len(i), i.index("*")+1+field)] for i in predict_list if i.startswith(gene)]))

    match_pair.sort(key=lambda i: i[1] or i[2])
    print(f"{'Answer':18} -- {'Predict':18}")
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


name = "linnil1_syn_wide"


def getSummaryAnswer(id):
    data = pd.read_csv(f"{name}/summary.csv", sep="\t", header=None)
    data = data[data[0].str.contains(id)]
    alleles = list(data[2].str.split("_"))[0]
    return alleles


def getPingResult(id):
    data = pd.read_csv(f"/home/linnil1/kir/PING/PING_test/{name}_result/finalAlleleCalls.csv", sep=',')
    data = data[data["Unnamed: 0"].str.contains(id)]
    data = data.drop(columns="Unnamed: 0").to_dict()
    alleles = []
    for i in data.values():
        for j in i.values():
            alleles.extend(j.split(' ')[0].split('+'))

    alleles = [i for i in alleles if "null" not in i]
    return alleles


def getHisatKIRResult(id, suffix="noab.sec_pair.noNH.tmp.group.nodup"):
    # suffix = "noab.sec_pair.noNH.tmp.group"
    data = open(f"data/{name}.{id}.{suffix}.txt").read().strip()
    alleles = data.split(',')[1:]
    return alleles



if __name__ == "__main__":
    # id = "01"
    data_result = []
    for i in range(0, 10):
        id = f"{i:02d}"
        print(id)
        data_result.append(evalute(getSummaryAnswer(id), getPingResult(id)))
        suffix="noabmsa.sec_pair.noNH.tmp.group"
        suffix="noab.sec_pair.noNH.tmp.group"
        suffix="noab.sec_pair.noNH.tmp.group.nodup"
        suffix="noabmsa.sec_pair.noNH.tmp.group.nodup"
        #data_result.append(evalute(getSummaryAnswer(id), getHisatKIRResult(id, suffix)))
    # pprint(data_result)

    TP, FP, FN, match_gene, match_3, total, cn_error = 0, 0, 0, 0, 0, 0, 0
    fail_allele, fail_sample = 0, 0
    for result in data_result:
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
