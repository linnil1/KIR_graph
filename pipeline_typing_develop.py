import numpy as np
import json
from collections import defaultdict
from pprint import pprint
# import matplotlib.pyplot as plt
import glob
import sys
import os
import re


def my_typing(allele_count_per_read):
    print("Start typing")

    top_n = 10
    s = defaultdict(int)
    groups = []

    if not len(allele_count_per_read):
        return [], []

    if len(allele_count_per_read) < 100:  # min reads
        return [], []

    id_map = set()
    for j in allele_count_per_read:
        id_map = id_map | dict(j['p']).keys() | dict(j['n']).keys()
    id_map = {i: ind for ind, i in enumerate(id_map)}

    allele_prob_per_read = []
    for j in allele_count_per_read:
        p, n = dict(j['p']), dict(j['n'])
        arr = np.ones(len(id_map))
        # s = {}
        # for i in p.keys() | n.keys():
        #     s[i] = 1.
        # for i in p:
        #     s[i] *= 0.99 ** p[i]
        # for i in n:
        #     s[i] *= 0.01 ** n[i]
        for i in p:
            arr[id_map[i]] *= 0.999 ** p[i]
        for i in n:
            arr[id_map[i]] *= 0.001 ** n[i]
        allele_prob_per_read.append(arr)

    allele_prob_per_read = np.array(allele_prob_per_read)

    # one allele
    allele_prob_all_1 = np.log(allele_prob_per_read).sum(axis=0) \
                        / len(allele_count_per_read)
    allele_prob_all_1_list = [(allele_prob_all_1[v1], [k1]) \
            for k1, v1 in id_map.items() if not k1.endswith("BACKBONE")]
    allele_prob_all_1_list = list(map(list, 
        sorted(allele_prob_all_1_list, key=lambda i:-i[0])))
    # pprint(allele_prob_all_1_list[:3])
    for max_allele in allele_prob_all_1_list[:top_n]:
        max_allele.append([np.mean(allele_prob_per_read.max(axis=1) == \
                                   allele_prob_per_read[:, id_map[max_allele[1][0]]]),
                           np.sum (allele_prob_per_read.max(axis=1) == \
                                   allele_prob_per_read[:, id_map[max_allele[1][0]]])])
        print(max_allele)
        
    # two allele
    allele_prob_all_2 = np.max([
        np.tile(allele_prob_per_read[:, None, :], [len(id_map), 1]),
        np.tile(allele_prob_per_read[:, :, None], [1, len(id_map)])], axis=0)
    allele_prob_all_2 = np.log(allele_prob_all_2).sum(axis=0) \
                        / len(allele_count_per_read)

    allele_prob_all_2_list = [(allele_prob_all_2[v1, v2], [k1, k2]) \
            for k1, v1 in id_map.items() \
                for k2, v2 in id_map.items() \
                    if (k1 <= k2 and
                        not k1.endswith("BACKBONE") and 
                        not k2.endswith("BACKBONE"))]
    allele_prob_all_2_list = list(map(list, 
        sorted(allele_prob_all_2_list, key=lambda i:-i[0])))

    for max_allele in allele_prob_all_2_list[:top_n]:
        allele_select = allele_prob_per_read[:, [id_map[max_allele[1][0]], id_map[max_allele[1][1]]]]
        nosame_ind = allele_select[:, 0] != allele_select[:, 1]
        allele_select_count = np.argmax(allele_select[nosame_ind], axis=1)
        max_allele.append([np.sum(allele_select_count == 0) + np.logical_not(nosame_ind).sum() // 2,
                           np.sum(allele_select_count == 1) + np.logical_not(nosame_ind).sum() // 2,
                           np.logical_not(nosame_ind).sum()])

        print(max_allele,
              sum(allele_select[:, 0] == allele_select[:, 1]))
    # pprint(allele_prob_all_2_list[:10])
    """
    return {'len': len(allele_count_per_read),
            '1':   allele_prob_all_1_list,
            '2':   allele_prob_all_2_list }
    """

    # third
    """
    first_n = min(100, len(allele_prob_all_2_list))
    allele_prob_all_3 = np.array([
        np.max([allele_prob_per_read[:, id_map[i[1]]],
                allele_prob_per_read[:, id_map[i[2]]]], axis=0)
            for i in allele_prob_all_2_list[:first_n]]).T
    allele_prob_all_3 = np.max([
        np.tile(allele_prob_per_read[:, None, :], [first_n,   1]),
        np.tile(allele_prob_all_3[:, :, None   ], [1, len(id_map)])], axis=0)
    allele_prob_all_3 = np.log(allele_prob_all_3).sum(axis=0) \
                        / len(allele_count_per_read)
    allele_prob_all_3_list = [(allele_prob_all_3[v1, v2], sorted([*allele_prob_all_2_list[v1][1], k3])) \
            for v1 in range(first_n) \
                for k3, v2 in id_map.items() \
                    if not k3.endswith("BACKBONE")]
    allele_prob_all_3_list = list(map(list, 
        sorted(allele_prob_all_3_list, key=lambda i:-i[0])))
    # pprint(allele_prob_all_3_list[:10])

    for max_allele in allele_prob_all_3_list[:3]:
        ids = [id_map[max_allele[i]] for i in range(1, 4)]
        allele_select = allele_prob_per_read[:, ids]
        v = (allele_select[:, 0] == allele_select[:, 1]) & (allele_select[:, 2] == allele_select[:, 1])
        allele_select_count = np.argmax(allele_select[np.logical_not(v)], axis=1)
        # max_allele.append(np.sum(allele_select_count == 0))
        # max_allele.append(np.sum(allele_select_count == 1))
        # max_allele.append(np.sum(allele_select_count == 2))
        print(max_allele,
              np.sum(allele_select_count == 0),
              np.sum(allele_select_count == 1),
              np.sum(allele_select_count == 2),
              sum(v))
    """
    # filter
    allele_prob_all_1_list = [i for i in allele_prob_all_1_list[:top_n] if i[2][0] > 0.9]
    allele_prob_all_2_list = [i for i in allele_prob_all_2_list[:top_n]
            if np.min(i[2][:2]) / np.max(i[2][:2]) > 0.3 and np.min(i[2][:2]) / i[2][2] > 0.2
                             ]

    # merge
    allele_prob_list = sorted([
                *allele_prob_all_1_list,
                *map(lambda i: [i[0] * 1.5, *i[1:]], allele_prob_all_2_list),
                # *map(lambda i: [i[0] * 1.5, i[1:]], allele_prob_all_3_list)
            ],
            key=lambda i: -i[0])

    pprint(allele_prob_list[:5], width=120)

    groups = []
    max_allele = allele_prob_list[0][1]
    allele_select = allele_prob_per_read[:, [*[id_map[a] for a in max_allele]]]
    groups = list(map(lambda i: max_allele[i], np.argmax(allele_select, axis=1)))

    return allele_prob_list, groups


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


# my_typing(json.load(open("data/GU182340_paired.pair.tmp.KIR2DL1.json")))
genes  = ["KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5",
          "KIR2DP1", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4",
          "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DP1",
          "KIR3DS1"]

# names = list(set(map(lambda i: i.split(".KIR")[0], glob.glob("data/synSeq.hiseq.dp50.rl150.*pair.tmp.*.json"))))
# print(names)
# exit()
name = f"data/synSeq.hiseq.dp50.rl150.1.split.pair.tmp"
# name = sys.argv[1]

name = f"data/synSeq.hiseq.dp50.rl150.1.merge.pair.tmp"
genes  = ["KIR"]

# reads all records in bam
headers = []
records_per_gene = defaultdict(list)
records_new = []
f_sam = open(f"{name}.sam")
for i in f_sam:
    if i[0] == "@":
        headers.append(i)
    else:
        gene = re.findall("\t(KIR.*?)\*BACKBONE\t", i)[0]
        records_per_gene[gene].append(i.strip())
        

# for gene in ["KIR3DL3"]:
for gene in genes:
    print(name, gene)
    data, groups = my_typing(json.load(open(name + f'.{gene}.json')))
    json.dump(data, open(name + f'.{gene}.em.json', 'w'), cls=MyEncoder)

    for i in set(groups):
        headers.append(f"@RG\tID:{i}\n")
    if len(groups) > 0:
        assert len(records_per_gene[gene]) == 2 * len(groups)
    for i, allele_name in enumerate(groups):
        records_new.append(records_per_gene[gene][i    ] + "\tRG:Z:" + allele_name + "\n")
        records_new.append(records_per_gene[gene][i + 1] + "\tRG:Z:" + allele_name + "\n")
    # os.system(f"samtools view -H {name}.sam > {name}.group.sam")

f_sam = open(f"{name}.group.sam", "w")
f_sam.writelines(headers)
f_sam.writelines(records_new)
f_sam.close()

os.system(f"samtools sort -@30 {name}.group.sam -o {name}.group.bam")
os.system(f"samtools index {name}.group.bam")



"""
pprint(Gene_counts[:100])
# [['KIR2DL1*002', 17496], ['KIR2DL1*061', 17004], ['KIR2DL1*062', 16948],
print(sorted(read_stat.items(), key=lambda i:-i[1])[:10])
[(0, 19302), (13, 2860), (2, 1873), (8, 1759), (24, 1566), (14, 1554), (11, 1488), (23, 1456), (10, 1337), (15, 1177)]
  (1, 632),

pprint(allele_stat)
print(sorted(allele_stat.items(), key=lambda i:-i[1])[:10])
 [(0, 9018), (13, 2781), (2, 1853), (8, 1759), (14, 1548), (24, 1543), (11, 1488), (23, 1456), (10, 1337), (15, 1167)]

# using allele
{'KIR2DL5B': 21186, 'KIR3DP1': 212305, 'KIR2DL1': 186795, 'KIR3DL1': 33365, 'KIR2DL2': 6117, 'KIR2DS2': 88234, 'KIR2DL3': 8438, 'KIR3DS1': 2430, 'KIR2DL4': 22062, 'KIR3DL2': 17634, 'KIR3DL3': 128934, 'KIR2DS4': 14478, 'KIR2DS3': 21746, 'KIR2DL5A': 16504, 'KIR2DP1': 16829, 'KIR2DS1': 23312, 'KIR2DS5': 27632}
[(1, 24036), (0, 9018), (2, 4442), (17, 1070), (3, 1061), (6, 25), (8, 11), (4, 10), (5, 8), (10, 6)]

# unique aligned allele read count
{'KIR2DL5B': 445, 'KIR2DL1': 5418, 'KIR3DL1': 456, 'KIR3DP1': 5328, 'KIR3DL2': 777, 'KIR2DS4': 559, 'KIR3DL3': 4205, 'KIR2DL3': 537, 'KIR2DS2': 3951, 'KIR2DP1': 1106, 'KIR2DS1': 2, 'KIR2DL2': 1242, 'KIR2DS5': 2, 'KIR2DS3': 1, 'KIR3DS1': 1, 'KIR2DL4': 6}

from Bio import SeqIO
alleles = SeqIO.to_dict(SeqIO.parse("kir_merge.save.fa", "fasta"))
kir_length = {}
for i in alleles:
    g = i.split("*")[0]
    if g not in kir_length:
        kir_length[g] = len(str(alleles[i].seq).replace("-", ""))

count = {'KIR2DL5B': 21186, 'KIR3DP1': 212305, 'KIR2DL1': 186795, 'KIR3DL1': 33365, 'KIR2DL2': 6117, 'KIR2DS2': 88234, 'KIR2DL3': 8438, 'KIR3DS1': 2430, 'KIR2DL4': 22062, 'KIR3DL2': 17634, 'KIR3DL3': 128934, 'KIR2DS4': 14478, 'KIR2DS3': 21746, 'KIR2DL5A': 16504, 'KIR2DP1': 16829, 'KIR2DS1': 23312, 'KIR2DS5': 27632}
for i in count:
    print(i, count[i] / kir_length[i] * 150)

# ensure pair-end
Reads: 39770 Pairs: 13081
{'KIR3DP1': 2537, 'KIR2DL1': 2668, 'KIR3DL1': 140, 'KIR2DS2': 1854, 'KIR2DL2': 540, 'KIR3DL3': 1546, 'KIR2DL3': 108, 'KIR2DS4': 88, 'KIR2DP1': 229, 'KIR3DL2': 381, 'KIR2DL5B': 242, 'KIR2DL4': 3, 'KIR2DS1': 3, 'KIR2DS5': 2}
[(0, 11680), (13, 2141), (8, 1192), (11, 1071), (2, 1054), (24, 1011), (15, 841), (23, 836), (14, 830), (10, 800)]
[(1, 10341), (2, 1134), (0, 986), (17, 335), (3, 259), (6, 10), (4, 4), (7, 3), (8, 3), (10, 3)]

# my method
defaultdict(<class 'int'>, {'KIR3DP1': 2890, 'KIR2DL1': 150, 'KIR2DP1': 428, 'KIR2DS1': 1601, 'KIR2DL2': 1766, 'KIR2DS2': 699, 'KIR3DL3': 1584, 'KIR3DS1': 63, 'KIR2DL3': 98, 'KIR2DS4': 177, 'KIR3DL2': 126, 'KIR3DL1': 38, 'KIR2DS3': 4, 'KIR2DS5': 10, 'KIR2DL4': 1})
[('ok', 9635), ('Cannot Determine', 2256), ('No allele found', 1190)]

# editdist > 5
Reads: 63129 Pairs: 27253
defaultdict(<class 'int'>, {'KIR3DP1': 4187, 'KIR2DL1': 883, 'KIR2DS4': 1013, 'KIR2DS2': 1930, 'KIR2DP1': 1483, 'KIR2DS1': 2541, 'KIR2DL2': 3194, 'KIR2DL3': 595, 'KIR3DL3': 3536, 'KIR3DS1': 158, 'KIR3DL2': 347, 'KIR3DL1': 149, 'KIR2DS3': 15, 'KIR2DL4': 2, 'KIR2DS5': 14})
[('ok', 20047), ('Cannot Determine', 4821), ('No allele found', 2385)]

# remove NH
Reads: 66749 Pairs: 28914
defaultdict(<class 'int'>, {'KIR3DP1': 4422, 'KIR2DL1': 900, 'KIR2DS4': 1088, 'KIR2DS2': 2003, 'KIR2DP1': 1547, 'KIR2DS1': 2708, 'KIR2DL2': 3409, 'KIR2DL3': 608, 'KIR3DL3': 3892, 'KIR3DS1': 169, 'KIR3DL2': 351, 'KIR3DL1': 164, 'KIR2DS3': 16, 'KIR2DL4': 2, 'KIR2DS5': 14})
[('ok', 21293), ('Cannot Determine', 5140), ('No allele found', 2481)]


def sort_group_order():
    os.system(f"samtools sort {alignment_fname} -o {bam_group_fname}")
    os.system(f"samtools index {bam_group_fname}")

    for gene, _ in Gene_counts:
        bam_group_split_name = bam_group_fname[:-4] + f".{gene}"
        os.system(f"samtools view -h -r {gene} {bam_group_fname} -o {bam_group_split_name}.bam")
        os.system(f"samtools index {bam_group_split_name}.bam")


# sort_group_order()
# allele_count_per_read = json.load(open("IHW01046.json"))
"""
