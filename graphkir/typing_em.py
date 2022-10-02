"""
The HISAT2 EM parts

* Read variants from JSON
* Allele calling by EM
"""
import sys
import json
from typing import TextIO
from itertools import chain
from collections import Counter, defaultdict
from dataclasses import dataclass

from Bio import SeqIO
import numpy as np
from .hisat2 import ReadsAndVariantsData, loadReadsAndVariantsData, removeMultipleMapped


@dataclass
class Hisat2AlleleResult:
    """ The result of hisat2 EM result """
    allele: str  # allele name
    count: int   # allele count
    prob: float  # allele abundance
    cn: int = 0  # Copy number (used in.kir_typing)


def readAlleleLength(file_fasta: str) -> dict[str, int]:
    """ Read allele's length """
    return {
        seq.id: len(seq.seq) for seq in SeqIO.parse(file_fasta, "fasta")
    }


def preprocessHisatReads(reads_data: ReadsAndVariantsData
                         ) -> dict[str, list[dict[str, list[list[str]]]]]:
    """
    Preprocess hisat2 reads before EM

    * Remove multiple mapped reads
    * Group by reference
    * map ID to variant
    * extract alleles from variant

    Returns:
        Dictionary[ backbone_name, list of postive and negative alleles of each read]
    """
    map_variant_to_alleles = {v.id: v.allele for v in reads_data['variants']}
    backbones_reads_alleles = defaultdict(list)
    pairreads = reads_data['reads']
    assert all(map(lambda i: i.multiple == 1, pairreads))

    for read in pairreads:
        backbones_reads_alleles[read.backbone].append({
            'lp': [map_variant_to_alleles[v] for v in read.lpv],
            'ln': [map_variant_to_alleles[v] for v in read.lnv],
            'rp': [map_variant_to_alleles[v] for v in read.rpv],
            'rn': [map_variant_to_alleles[v] for v in read.rnv],
        })
    return backbones_reads_alleles


def getCandidateAllelePerRead(positive_allele: list[list[str]],
                              negative_allele: list[list[str]]) -> list[str]:
    """
    positive_allele & positive_allele - negative_allele - negative_allele

    same as get_count in hisatgenotype
    """
    candidate = None
    for allele in positive_allele:
        if candidate is None:
            candidate = set(allele)
        else:
            candidate &= set(allele)
    if candidate is None:
        return []

    for allele in negative_allele:
        candidate -= set(allele)
    return list(candidate)


def getMostFreqAllele(candidates: list[str]) -> list[str]:
    """
    Find the most frequent alleles in those variants

    It may return multiple alleles that all of them are in high frequency

    same as get_stat in hisatgenotype
    """
    # count the allele
    count = Counter(candidates)
    if not count:
        return []
    max_count = max(count.values())
    max_allele = [i[0] for i in count.items() if i[1] == max_count]
    return max_allele


def hisatEMnp(
        allele_per_read: list[list[str]],
        seq_len: dict[str, int] = {},
        iter_max: int = 300,
        diff_threshold: float = 0.0001,
        ) -> dict[str, float]:
    """
    EM-algorithm: allele per reads -> allele probility

    same as single_abundance in hisat-genotype

    Accelerated version of EM - SQUAREM iteration
    Varadhan, R. & Roland, C. Scand. J. Stat. 35, 335-353 (2008)
    Also, this algorithm is used in Sailfish
    http://www.nature.com/nbt/journal/v32/n5/full/nbt.2862.html


    Args:
        allele_per_read: The list of alleles per reads
        seq_len:
             The length of the sequences
             (If not provided, it'll not normalized by length)
        iter_max: Max iteration
        diff_threshold: The threshold to stop iteration (The difference between each step)

    Returns:
        The abundance of each allele
    """
    # init probility (all allele has even probility)
    allele_name = sorted(set(chain.from_iterable(allele_per_read)))
    allele_map = dict(zip(allele_name, range(len(allele_name))))

    if seq_len:
        allele_len = np.array([seq_len[i] for i in allele_name])
    else:
        allele_len = np.ones(len(allele_name))
    allele_select_one = []
    for alleles in allele_per_read:
        x = np.zeros(len(allele_map))
        for i in map(allele_map.get, alleles):
            x[i] = 1
        allele_select_one.append(x)
    allele_select_one = np.array(allele_select_one)  # type: ignore

    def getNextProb(prob_per_allele):
        a = prob_per_allele * allele_select_one
        b = a.sum(axis=1)[:, None]
        a = np.divide(a, b,
                      out=np.zeros(a.shape),
                      where=b != 0)
        a /= allele_len
        a = a.sum(axis=0)
        return a / a.sum()

    # Run EM
    prob = getNextProb(np.ones(len(allele_map)))
    for iters in range(0, iter_max):
        # SQUAREM
        prob_next = getNextProb(prob)
        prob_next2 = getNextProb(prob_next)
        r = prob_next - prob
        v = prob_next2 - prob_next - r
        r_sum = (r ** 2).sum()
        v_sum = (v ** 2).sum()
        if v_sum > 0.0:
            g = -np.sqrt(r_sum / v_sum)
            prob_next3 = prob - r * g * 2 + v * g ** 2
            prob_next3 = np.maximum(prob_next3, 0)
            prob_next = getNextProb(prob_next3)

        # Recalculate diff between previous stage
        diff = np.abs(prob - prob_next).sum()
        if diff <= diff_threshold:
            break

        # next
        # print(f"Iter: {iters} Diff: {diff}")
        prob = prob_next

    return dict(zip(allele_name, prob))


def hisat2TypingPerGene(reads_alleles: list[dict[str, list[list[str]]]]
                        ) -> list[Hisat2AlleleResult]:
    """  The orignal typing method in hisat2 (pre gene) """
    reads_max_alleles = []
    for read in reads_alleles:
        reads_max_alleles.append(
            getMostFreqAllele(
                getCandidateAllelePerRead(read['lp'], read['ln']) +
                getCandidateAllelePerRead(read['rp'], read['rn'])
            )
        )

    allele_prob = hisatEMnp(reads_max_alleles)
    allele_count = Counter(chain.from_iterable(reads_max_alleles))

    alleles_stat = [Hisat2AlleleResult(
        allele=allele,
        count=allele_count[allele],
        prob=allele_prob[allele],
    ) for allele in allele_prob.keys() | allele_count.keys()]
    return alleles_stat


def hisat2Typing(read_and_variant_json: str, output_prefix: str):
    """
    The orignal typing method in hisat2

    Args:
      read_and_variant_json: The json file in `gk_hist2.py`
      output_prefix: The result is save in json and text (report)
    """
    reads_data = loadReadsAndVariantsData(read_and_variant_json)
    reads_data = removeMultipleMapped(reads_data)
    hisat_result = {}

    for backbone, reads_alleles in preprocessHisatReads(reads_data).items():
        hisat_result[backbone] = hisat2TypingPerGene(reads_alleles)

    printHisatTyping(hisat_result)
    with open(output_prefix + ".txt", "w") as f:
        printHisatTyping(hisat_result, file=f)
    with open(output_prefix + ".json", "w") as f:
        json.dump(hisat_result, f)


def printHisatTyping(hisat_result: dict[str, list[Hisat2AlleleResult]],
                     first_n: int = 10,
                     file: TextIO = sys.stdout):
    """
    Print the typing result (EM)

    Args:
        hisat_result: Hisat2AlleleResult grouped by backbone
        first_n: Print top-n alleles
        file: Output text to the handle (default: stdout)
    """
    for backbone, result in hisat_result.items():
        print(backbone, file=file)
        allele_count = sorted(result, key=lambda i: i.count, reverse=True)
        for i, allele in enumerate(allele_count[:first_n]):
            print(f"  {i+1:2d} {allele.allele:18s} "
                  f"(count: {allele.count})", file=file)

        allele_prob = sorted(result, key=lambda i: i.prob, reverse=True)
        for i, allele in enumerate(allele_prob[:first_n]):
            print(f"  Rank {i+1:2d} {allele.allele:18s} "
                  f"(abundance: {allele.prob:.2f}, cn: {allele.cn})", file=file)
