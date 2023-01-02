# This script is inspired from hisat-genotype/hisatgenotype_typing_core
import re
import os
import sys
import json
import copy
import bisect
from typing import Union, ClassVar
from collections import defaultdict, Counter
from dataclasses import dataclass, field
import numpy as np
from Bio import SeqIO

from kg_utils import runDocker, getSamples

# bubble
import matplotlib
matplotlib.use("webagg")
matplotlib.rcParams['webagg.port'] = 8082
matplotlib.rcParams['webagg.open_in_browser'] = False
import matplotlib.pyplot as plt
import networkx as nx
import pickle

# sam
num_editdist = 4
# sam
NH_filtering = False
# em
diff_threshold = 0.0001
# em
iter_max = 100
# em
norm_by_length = False
# pileup
error_correction = False
# output
extend_allele = True # set false if merge(too big)

from kg_typping import readPair, filterPair, record2Variant, getVariantsBoundary
import kg_typping


no_novel = False
kg_typping.error_correction = error_correction


@dataclass
class Variant:
    pos: int
    typ: str
    ref: str = None
    val: Union[int, str] = None
    allele: list = field(default_factory=list)
    length: int = 1
    in_exon: bool = False
    id: str = None

    # allele: list = field(default_factory=list)
    # freq: float = None
    order_type: ClassVar[dict] = {"insertion": 0, "single": 1, "deletion": 2}
    order_nuc: ClassVar[dict] = {"A": 0, "C": 1, "G": 2, "T": 3}
    novel_count: ClassVar[int] = 0

    def __lt__(self, othr):
        # The sorted sta
        return (self.ref, self.pos) < (othr.ref, othr.pos)

    def __eq__(self, othr):
        return (self.pos, self.ref, self.typ, self.val) \
               == (othr.pos, othr.ref, othr.typ, othr.val)

    def __hash__(self):
        return hash((self.pos, self.ref, self.typ, self.val))

    def __repr__(self):
        return self.id


def findVariantId(variant):
    global novel_variants, variants
    """ Match variant to hisat2 .snp """
    # The variant is not located in .snp.index but in .snp
    if variant in variants:
        # print("Find it", variant)
        return variants[variant]
    elif variant in novel_variants:
        return novel_variants[variant]
    # cannot find: set it novel
    elif variant.typ in ["single", "insertion", "deletion"]:
        # print("Cannot find", variant)
        variant.id = f"nv{Variant.novel_count}"
        Variant.novel_count += 1
        novel_variants[variant] = variant
        return variant
    else:
        # print("Other", variant)
        assert variant.typ == "match"
        return variant


def getPNFromVariantList(variant_list, exon_only=False):
    """
    Extract **useful** positive and negative varinat from the list

    Note:
      * Remove novel insertion and deletion
      * Find negative variant (exclude deletion at the front and end)
    """
    # Find all variant in the region
    variant_list = list(variant_list)
    if not variant_list:
        return [], []
    left, right = getVariantsBoundary(variant_list)
    pos_right = variant_list[-1].pos + variant_list[-1].length
    assert left <= right

    # TODO: linnil1
    # insert or novel deletion = mapping error
    # if any(v.typ == "insertion" for v in variant_list):
    #     return [], []
    # if any(v.typ == "deletion" and v.id.startswith("nv") for v in variant_list):
    #     return [], []

    # exclude for negative
    exclude_variants = set()

    # exclude cannot error correction variant
    for v in variant_list:
        if v.val == "N":
            for i in "ATCG":
                va = copy.deepcopy(v)
                va.val = i
                exclude_variants.add(va)

    # remove match
    variant_list = [v for v in variant_list if v.typ != "match"]
    # remove novel
    if no_novel:
        variant_list = [v for v in variant_list if v.id and not v.id.startswith("nv")]

    # positive variation
    if exon_only:
        positive_variants = [v for v in variant_list if v.in_exon]
    else:
        positive_variants = variant_list
    exclude_variants.update(positive_variants)

    # print('positive', variant_list)

    # negative
    negative_variants = []
    for v in variants_sorted_list[left:right]:
        if v in exclude_variants:
            continue
        if exon_only and not v.in_exon:
            continue
        # Deletion should not over the right end of read
        # Because some deletion is ambiguous until last one
        # TODO: change 10 -> to repeat length
        # if v.typ == "deletion" and v.pos + v.val + 10 >= pos_right:
        #     continue
        # print('negative', v)
        negative_variants.append(v)
        # TODO: maybe some deletion is before left
        # they use gene_var_maxrights to deal with

    return positive_variants, negative_variants


def recordToVariants(record, pileup):
    variant_list, soft_clip = record2Variant(record)
    # linnil1: soft clip = mapping error
    # if sum(soft_clip) > 0:
    #     return []
    # variant_list = flatten(map(lambda i: pileupModify(pileup, i), variant_list))
    variant_list = map(findVariantId, variant_list)
    # no sure what's this doing
    # variant_list = extendMatch(variant_list)  # remove novel -> extend
    return variant_list


# index
def setIndex(index):
    global novel_variants
    novel_variants = {}
    global variants, variants_sorted_list
    variants = list(kg_typping.variants.keys())
    # bubble tmp
    max_pos = max(v.pos for v in variants)
    refs = set(v.ref for v in variants)
    for i in range(10, max_pos, 100):
        for r in refs:
            variants.append(Variant(typ="ref",
                                    ref=r,
                                    pos=i,
                                    length=1,
                                    id=f"p{i // 100}",
                                    val=r))
    variants = {v: v for v in variants}
    variants_sorted_list = sorted([i for i in variants])
    kg_typping.variants = variants
    kg_typping.variants_sorted_list = variants_sorted_list


def flatten(its):
    for it in its:
        yield from it


def pileupModify(pileup, variant):
    if variant.typ == "single":
        p = pileup.get((variant.ref, variant.pos))
        if not p:
            # print("error pileup", variant)
            yield variant
            return

        # error correction critria
        if p['all'] > 20:
            if p.get(variant.val, 0) > 0.2:
                # print("reserve", variant, p)
                yield variant
                return

            # case: AAAAAAAATAAAA where REF = G
            if any(i[1] >= 0.8 for i in p.items() if i[0] != "all"):
                variant.val = max([i for i in p.items() if i[0] != "all"], key=lambda i: i[1])[0]
                yield variant
                return

            # if variant.id.startswith("hv"):
            #     print("to_match", variant, p)
            # case: AAAAAACCCCCCCCCCC REF = C
            # neither set to A or C is wrong
            # if I set it to match or novel -> it will become negative
            # variant.typ = "match"
            variant.val = "N"
            yield variant

        # else -> insufficient info
        yield variant
        return

    # TODO: split the match when some base in match is sequencing error
    yield variant
    return


def main(bam_file):
    new_suffix = ".bubble"
    # if exon_only: new_suffix += ".exon"
    if error_correction:
        # pileup
        new_suffix += ".errcorr"
        print("Reading pileup")
        pileup = getPileupDict(bam_file)
        print("Reading pileup Done")
    else:
        pileup = {}

    name_out = os.path.splitext(bam_file)[0] + new_suffix

    # TODO: test
    global extend_allele
    if "_merge" in bam_file:
        extend_allele = False

    # read bam file
    pair_reads = readPair(bam_file)
    pair_reads = filter(filterPair, pair_reads)

    save_reads = defaultdict(list)
    hisat_gene_alleles = defaultdict(list)
    tot_pair = 0

    plot_ref = "KIR2DL1S1*BACKBONE"
    # = ref -> level0
    # = alt -> level -1
    # = novel -> level 1
    G = nx.MultiGraph()
    pos_count = defaultdict(int)
    for i, v in enumerate(sorted(variants)):
        if v.ref != plot_ref:
            continue
        if v.typ != "ref":
            # negative variant
            rv = copy.deepcopy(v)
            rv.val = v.ref
            rv.id = "r" + v.id
            G.add_node(rv, pos=(v.pos + pos_count[v.pos] * 0.2, pos_count[v.pos]))
            pos_count[v.pos] -= 1

        # positive variant
        G.add_node(v, pos=(v.pos + pos_count[v.pos] * 0.2, pos_count[v.pos]))
        pos_count[v.pos] -= 1
    pos_count = defaultdict(lambda: 1)

    # same
    save_alleles = []

    # main
    for left_record, right_record in pair_reads:
        if f"\t{plot_ref}\t" not in left_record:
            continue
        tot_pair += 1
        # if tot_pair > 100:
        #     break

        # group read by gene name
        backbone = left_record.split("\t")[2]

        # if "KIR3DL3*0090101-1-2942" not in left_record:
        #     continue

        lv = recordToVariants(left_record, pileup)
        rv = recordToVariants(right_record, pileup)
        lpv, lnv = getPNFromVariantList(lv)
        rpv, rnv = getPNFromVariantList(rv)

        # negative -> ref
        lnv_ref = [copy.deepcopy(i) for i in lnv]
        for i in lnv_ref:
            i.val = i.ref
        rnv_ref = [copy.deepcopy(i) for i in rnv]
        for i in rnv_ref:
            i.val = i.ref

        # merge pos and negative variant
        l_v = sorted([*lpv, *lnv_ref])
        r_v = sorted([*rpv, *rnv_ref])

        if set(l_v) & set(r_v):
            print(l_v, r_v, set(l_v) & set(r_v))

        # set novel position
        for v in [*l_v, *r_v]:
            if v.id.startswith("nv"):
                G.add_node(v, pos=(v.pos + pos_count[v.pos] * 0.2, pos_count[v.pos]))
                pos_count[v.pos] += 1

        # add to graph
        for i in range(1, len(l_v)):
            v = l_v[i]
            G.add_edge(l_v[i], l_v[i-1])
        for i in range(1, len(r_v)):
            v = r_v[i]
            G.add_edge(r_v[i], r_v[i-1])

        save_alleles.append({
            'l_sam': left_record,
            'r_sam': right_record,
            'left': l_v,
            'right': r_v,
        })

    # graph -> weighted graph
    Gw = nx.Graph(G)
    for a, b, _ in G.edges(data=True):
        if 'weight' in Gw[a][b]:
            Gw[a][b]['weight'] += 1
        else:
            Gw[a][b]['weight'] = 1

    # filter small edge and 0 node
    def filterEdge(a, b):
        return Gw[a][b]['weight'] > 1
    Gw_filtered = Gw
    Gw_filtered = nx.Graph(nx.subgraph_view(Gw_filtered, filter_edge=filterEdge))
    Gw_filtered = nx.Graph(nx.subgraph_view(Gw_filtered, filter_node=lambda i: Gw_filtered.degree(i, weight="weight") > 2))

    # save
    name = "tmp"
    if no_novel:
        name += ".no_novel"
    nx.write_gpickle(Gw_filtered, f"{name}.pkl")
    pickle.dump(save_alleles, open(f"{name}.v.pkl", "wb"))

    # plot
    pos = nx.get_node_attributes(G, 'pos')
    nx.draw(Gw_filtered, pos, arrows=True, with_labels=True, font_size=8, font_color="red", edgelist=[])
    for a, b, weight in Gw_filtered.edges(data='weight'):
        print(a, b, weight)
        nx.draw_networkx_edges(Gw_filtered, pos, edgelist=[(a, b)], width=np.sqrt(weight), connectionstyle='arc3, rad = 0.3')
    plt.show()

    print("Filterd pairs", tot_pair)
    return new_suffix


if __name__ == "__main__":
    kg_typping.setIndex("index/kir_2100_2dl1s1.mut01")
    setIndex("index/kir_2100_2dl1s1.mut01")
    name = "data/linnil1_syn_30x.00.index_kir_2100_2dl1s1.mut01.bam"
    main(name)
