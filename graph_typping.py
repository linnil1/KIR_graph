import pickle
from collections import Counter
import numpy as np
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt

from kg_fine_typping import Variant
from kg_linnil1 import onehot2Prob, variant2onehot, calcProb, printProb
from pprint import pprint


# matplotlib.rcParams['backend'] = "QtCairo"
matplotlib.use("webagg")
matplotlib.rcParams['webagg.port'] = 8082
matplotlib.rcParams['webagg.open_in_browser'] = False
G = nx.read_gpickle("tmp.no_novel.pkl")
save_alleles = pickle.load(open("tmp.no_novel.v.pkl", 'rb'))
for allele in save_alleles:
    allele['left'] = set(allele['left'])
    allele['right'] = set(allele['right'])
# print(save_alleles)

# print(G.nodes(data=True))
# print(G.edges(data=True))


def flatten(list_list_data):
    for list_data in list_list_data:
        yield from list_data


def plotEdges(G):
    for a, b, weight in G.edges(data='weight'):
        # print(a,b)
        nx.draw_networkx_edges(G, pos, edgelist=[(a, b)], width=np.sqrt(weight), connectionstyle='arc3, rad = 0.3')


def getAlleles(variant):
    if variant.typ == 'ref':
        return {}
    # rhv xx
    if variant.ref == variant.val:
        return all_alleles - set(variant.allele)

    # hvxx
    return set(variant.allele)


nodes = G.nodes()
all_alleles = set(flatten(map(lambda i: i.allele, nodes)))
all_alleles_map = {name: i for i, name in enumerate(sorted(all_alleles))}
all_alleles_map_rev = {v: k for k, v in all_alleles_map.items()}
non_cycle_nodes = set(nodes) - set(flatten(nx.cycle_basis(G)))
# print(non_cycle_nodes)
# c = Counter(flatten(map(getAlleles, non_cycle_nodes)))
# print(c)
Gc = nx.subgraph_view(G, filter_node=lambda i: i not in non_cycle_nodes)
coms = list(nx.connected_components(Gc))


def calcReadAlleleProb(reads):
    reads_arr = []
    for variants in reads:
        read_arr = [onehot2Prob(variant2onehot(all_alleles_map, getAlleles(i))) for i in variants]
        reads_arr.append(np.stack(read_arr).prod(axis=0))
    reads = np.stack(reads_arr)
    return reads


def getAlleleSet():
    pass


def calcAlleleProb(reads):
    reads = calcReadAlleleProb(reads)
    log_probs = np.log(reads)

    # find the maximum probility of allele across all reads
    top_n = 100
    prob_save = []
    iter_max = 3
    for prob_iter in range(1, 1 + iter_max):
        if prob_save:
            prob_last = prob_save[-1]
        else:
            prob_last = None

        prob_res = calcProb(log_probs, prob_last,
                            prob_iter=prob_iter,
                            allele_length=np.ones(len(all_alleles_map)),
                            top_n=top_n)
        prob_res['allele_name'] = [list(map(all_alleles_map_rev.get, i)) for i in prob_res['allele_id']]
        prob_save.append(prob_res)
    # printProb(prob_save[-1])
    return [tuple(i) for i in prob_save[-1]['allele_name']]


i = 0
probs = []
for com in coms:
    # com_count = Counter(flatten(map(getAlleles, com)))
    reads = []
    for allele in save_alleles:
        if com & allele['left']:
            reads.append(com & allele['left'])
        if com & allele['right']:
            reads.append(com & allele['right'])
    # print(com, reads)
    probs.append(calcAlleleProb(reads))
    i += 1
    # if i > 10:
    #     break

print(i)
count = Counter(flatten(probs))
pprint(sorted(count.items(), key=lambda i: -i[1])[:100])
# for prob in probs:
#    printProb(prob)
exit()


# plot
pos = nx.get_node_attributes(G, 'pos')
# nx.draw(G, pos, nodelist=non_cycle_nodes, node_color="blue", edgelist=[])
# nx.draw(G, pos, arrows=True, with_labels=True, font_size=8, font_color="red", edgelist=[])

# with label will output all labels (not only in node_list)
nx.draw(G, pos, nodelist=non_cycle_nodes, node_color="blue", font_size=8, font_color="red", with_labels=True, edgelist=[])
for com in coms:
    print(com)
    nx.draw_networkx_nodes(Gc, pos, nodelist=com)
plotEdges(Gc)
plt.show()
