"""
Our main typing method
"""
from __future__ import annotations
import copy
from typing import Optional
from itertools import chain
from collections import defaultdict
from dataclasses import dataclass, field
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from .hisat2 import PairRead, Variant


@dataclass
class TypingResult:
    """
    Our typing result of each CN step

    Attributes:
      n: CN
      value:       The log-likelihood of the top-n allele set
      allele_id:   The id in top-n allele set
      allele_name: The allele names in top-n allele set
      allele_prob: The probility of the top-n allele set for each read
      fraction:    The perpotion of reads belonging to the allele in allele set
      value_sum_indv: The value of each allele at CN=1
      allele_name_group:
          The group of alleles that belongs to each allele in allele_name(s)
          e.g.
          ```
          allele_name = [00201, 00101]
          allele_name_group = [[0020101, 0020102], [00106, 0010103]]
          ```
    """
    n: int
    value: np.ndarray             # size: top_n
    value_sum_indv: np.ndarray    # size: top_n x n
    allele_id: np.ndarray         # size: top_n x n
    allele_name: list[list[str]]  # size: top_n x n
    allele_prob: np.ndarray       # size: read x top_n
    fraction: np.ndarray          # size: top_n x n
    allele_name_group: list[list[list[str]]] = field(default_factory=list)
                                  # size: top_n x n x group_size

    def selectBest(self) -> list[str]:
        """
        Select the best allele set by the maximum of likelihood
        whlie considering abundance (fraction).

        The criteria is that the allele set will be ignore
        if one the allele has abundance small than 0.5 / CN.

        If all allele set cannot meet the criteria, return first result.

        Example:
          loss -1  fraction 0.1 0.9
          loss -2  fraction 0.05 0.95
          loss -3  fraction 0.2 0.8
          loss -4  fraction 0.4 0.6  -> OK
        """
        best_id = 0
        for i in range(len(self.fraction)):
            expect_prob = 1 / self.n
            if all(j >= expect_prob / 2 for j in self.fraction[i]):
                break
            best_id += 1
        else:  # cannot meet the criteria
            best_id = 0
        return self.allele_name[best_id]

    def print(self, num: int = 10):
        """
        Print the first `num` typing result

        Example:
            ```
            2
            0 -52.37587690948907
              id   1 name KIR3DS1*0130108      fraction 0.35294
              id   2 name KIR3DS1*078          fraction 0.64705
            1 -52.37587690948907
              id   0 name KIR3DS1*0130107      fraction 0.42156
              id   2 name KIR3DS1*078          fraction 0.57843
            ```
        """
        print(self.n)
        top_n = len(self.value)
        n = self.allele_id.shape[1]

        for rank in range(min(top_n, num)):
            print("Rank",      rank,
                  "probility", self.value[rank],
                  "sum",       self.value_sum_indv[rank].sum())

            for i in range(n):
                print("  id",       f"{self.allele_id[rank][i]:3}",       end=" ")
                print("  name",     f"{self.allele_name[rank][i]:20s}",   end=" ")
                print("  fraction", f"{self.fraction[rank][i]:.5f}",      end=" ")
                print("  sum",      f"{self.value_sum_indv[rank][i]:3f}", end=" ")
                # print(f" unique_count {unique_count}")
                if self.allele_name_group:
                    print("  group", f"{self.allele_name_group[rank][i]}", end=" ")
                print()

    def fillNameGroup(self, allele_group_mapping: dict[str, list[str]]):
        """ extend the allele groups by allele_name and save in allele_name_group """
        self.allele_name_group = [[allele_group_mapping[j] for j in i]
                                  for i in self.allele_name]

    def sortByScoreAndEveness(self) -> TypingResult:
        """ reorder the internal top_n axis by score """
        rank_index = rankScore(self.value, self.value_sum_indv, self.fraction)
        return TypingResult(
            n              = self.n,
            value          = self.value          [rank_index],
            value_sum_indv = self.value_sum_indv [rank_index],
            allele_id      = self.allele_id      [rank_index],
            allele_name    = [self.allele_name[i] for i in rank_index],
            allele_prob    = self.allele_prob [:, rank_index],
            fraction       = self.fraction       [rank_index],
        )


def argSortRow(data: np.ndarray) -> list[int]:
    """ Argsort the data by row """
    return sorted(range(len(data)), key=lambda i: tuple(data[i]))


def rankScore(value: np.ndarray, value_sum_indv: np.ndarray, fraction: np.ndarray) -> list[int]:
    """
    Sort by likelihood score, sum of likelihood score per alele,
    and difference of allele abundance (smaller -> evenly-distributed)
    """
    fraction_diff = fraction - fraction.mean(axis=1, keepdims=True)           # size: top_n x CN
    fraction_diff = np.abs(fraction_diff).sum(axis=1)                         # size: top_n
    rank_index    = argSortRow(np.array([-value,
                                         -value_sum_indv.sum(axis=1),
                                         fraction_diff]).T)
    return rank_index                                                         # size: top_n


class AlleleTyping:
    """
    Our proposed method: Allele typing for multiple alleles

    Attributes:
      top_n (int): Consider top-n maximum likelihood to reduce computation
      results (list[TypingResult]): The typing result per steps
      probs (np.ndarray): The probility of read belong to the allele/allele-set
      log_probs (np.ndarray): log10 of probs
      id_to_allele (dict[int, str]): map the allele id to allele name
    """

    def __init__(self, reads: list[PairRead], variants: list[Variant]):
        """
        Parameters:
          reads: The reads belonged to the gene
          variants:
            The variants list.
            The variant_id in each reads will map to the variant for
            getting related alleles.
        """
        self.top_n = 300

        self.variants: dict[str, Variant] = {str(v.id): v for v in variants}
        allele_names = self.collectAlleleNames(variants)
        self.id_to_allele: dict[int, str] = dict(enumerate(sorted(allele_names)))
        self.allele_to_id: dict[str, int] = {j: i for i, j in self.id_to_allele.items()}

        self.probs = self.reads2AlleleProb(self.removeEmptyReads(reads))
        self.log_probs = np.log10(self.probs)

        self.result: list[TypingResult] = []
        # , allele_length: dict[str, float]):
        # allele_length = np.array([allele_length_map[self.allele_name_map_rev[i]]
        # for i in range(len(self.allele_name_map_rev))])
        # self.allele_length = allele_length / 10000.  # just a non-important normalize

    @staticmethod
    def removeEmptyReads(reads: list[PairRead]) -> list[PairRead]:
        """
        If the reads don't contains any information,
        i.e. no positive/negative variants in the read,
        discard it
        """
        return [read for read in reads if read.lpv + read.lnv + read.rpv + read.rnv]

    @staticmethod
    def collectAlleleNames(variants: list[Variant]) -> set[str]:
        return set(chain.from_iterable(map(lambda i: i.allele, variants)))

    def read2Onehot(self, variant: str) -> np.ndarray:
        """ Convert allele names in the variant into onehot encoding """
        onehot = np.zeros(len(self.allele_to_id), dtype=bool)
        for allele in self.variants[variant].allele:
            onehot[self.allele_to_id[allele]] = True
        return onehot

    @staticmethod
    def onehot2Prob(onehot: np.ndarray) -> np.ndarray:
        """ Onehot encoding -> probility"""
        # TODO: use quality
        prob = np.ones(onehot.shape) * 0.001
        prob[onehot] = 0.999
        return prob

    def reads2AlleleProb(self, reads: list[PairRead]) -> np.ndarray:
        """ Position/Negative variants in read -> probility of read belonged to allele """
        probs = []
        for read in reads:
            prob = [
                *[self.onehot2Prob(               self.read2Onehot(i) ) for i in read.lpv],
                *[self.onehot2Prob(               self.read2Onehot(i) ) for i in read.rpv],
                *[self.onehot2Prob(np.logical_not(self.read2Onehot(i))) for i in read.lnv],
                *[self.onehot2Prob(np.logical_not(self.read2Onehot(i))) for i in read.rnv],
            ]
            if not prob:
                continue
            probs.append(np.stack(prob).prod(axis=0))
        # probs = [i for i in probs if i is not None]
        return np.stack(probs)

    def typing(self, cn: int) -> TypingResult:
        """
        Typing cn alleles.

        Returns:
          The top-n allele-set are best fit the reads (with maximum probility).
            Each set has CN alleles.
        """
        self.result = []
        for _ in range(cn):
            res = self.addCandidate()
            res.print()
        return self.result[-1]

    def mapAlleleIDs(self, list_ids: np.ndarray) -> list[list[str]]:
        """ id (m x n np array) -> name (m list x n list of str)"""
        return [[self.id_to_allele[id] for id in ids] for ids in list_ids]

    @staticmethod
    def uniqueAllele(data: np.ndarray) -> np.ndarray:
        """
        Unique the allele set, return the mask.

        The boolean in the mask: 1 for unqiue allele set 0 for non-unique.

        Example:
          * data [[0,1], [1,0], [2,0], [2, 2], [2,0], [1, 0]]
          * return [1, 0, 1, 1, 0, 0]
        """
        s = set()
        index_unique = []
        for ids in data:
            ids = tuple(sorted(ids))
            if ids in s:
                index_unique.append(False)
            else:
                s.add(ids)
                index_unique.append(True)
        return np.array(index_unique)

    def addCandidate(self, candidate_allele: Optional[list[str]] = None) -> TypingResult:
        """
        The step of finding the maximum likelihood allele-set that can fit all reads

        Once call this function, the number alleles in the allele set increate 1.
        And the result will be saved in `.result`.

        Parameters:
            candidate_allele: The candidate alleles for this iteration
                              set it None for selecting all alleles
        Returns:
            The result of this iteration
        """
        if candidate_allele is None:
            allele_index = np.arange(self.log_probs.shape[1])                  # size: allele (select all)
        else:
            allele_index = np.array(list(map(self.allele_to_id.get, candidate_allele)))
                                                                               # size: allele (select partial)
        if not len(self.result):
            # first time: i.e. CN = 1
            allele_1_prob      = self.log_probs[:, allele_index].sum(axis=0)   # size: allele (selected array)
            allele_1_id        = allele_index[:, None]                         # size: top_n x 1
            # sort by decending order
            allele_1_top_index = np.argsort(allele_1_prob)[::-1][:self.top_n]  # size: top_n (selected array index)
            allele_1_top_id    = allele_1_id[allele_1_top_index]               # size: top_n
            allele_1_top_prob  = self.log_probs[:, allele_1_top_id.flatten()]  # size: reads x top_n
            allele_1_top_value = allele_1_prob[allele_1_top_index]             # size: top_n

            self.result.append(TypingResult(
                n              = 1,
                value          = allele_1_top_value,
                value_sum_indv = allele_1_top_value[:, None],                  # size: top_n x 1
                allele_id      = allele_1_top_id,                              # size: top_n x 1
                allele_name    = self.mapAlleleIDs(allele_1_top_id),           # size: top_n x 1
                allele_prob    = allele_1_top_prob,
                fraction       = np.ones(allele_1_top_id.shape),               # size:  top_n x 1
            ))
            return self.result[-1]

        # get previous top-n allele
        allele_n_1_prob = self.result[-1].allele_prob  # size: read x top_n
        allele_n_1_id   = self.result[-1].allele_id    # size: top_n x CN

        # Maximum
        # (top-n x (allele-1 ... allele-n-1)) x (allele_1 ... allele_m)
        allele_n_prob = np.maximum(self.log_probs[:, allele_index],
                                   allele_n_1_prob.T[:, :, None])             # size: top_n x read x allele
        allele_n_prob = allele_n_prob.sum(axis=1).flatten()                   # size: (top_n * allele)

        # Mean
        # cn = len(self.result)
        # allele_n_prob = self.log_probs[:, allele_index] \
        #                 + allele_n_1_prob.T[:, :, None] * cn                  # size: top_n x read x allele
        # allele_n_prob = allele_n_prob.sum(axis=1).flatten() / (cn + 1)        # size: top_n

        # create id (size: top_n * allele x n) of the allele_n_prob
        allele_n_id = np.hstack([
            # [[1],[3],[5]] -> [[1],[3],[5],[1],[3],[5]]
            np.repeat(allele_n_1_id, len(allele_index), axis=0),
            # [1,3,5] -> [1,3,5,1,3,5]
            np.tile(allele_index, len(allele_n_1_id))[:, None],
        ])

        # unique the allele set
        unique_id_index    = self.uniqueAllele(allele_n_id)
        allele_n_id        = allele_n_id[unique_id_index]
        allele_n_prob      = allele_n_prob[unique_id_index]

        # top-n result
        allele_n_top_index = np.argsort(allele_n_prob)[::-1][:self.top_n]      # size: top_n
        allele_n_top_id    = allele_n_id[allele_n_top_index]                   # size: top_n x CN
        allele_n_top_prob  = self.log_probs[:, allele_n_top_id].max(axis=2)    # read x top_n x CN -> size: read x top-n
        allele_n_top_value = allele_n_prob[allele_n_top_index]                 # size: top_n
        allele_n_top_sum   = self.log_probs[:, allele_n_top_id].sum(axis=0)    # size: top_n

        # calculate fraction
        belong             = np.equal(self.log_probs[:, allele_n_top_id],
                                      allele_n_top_prob[:, :, None])           # size: reads x top-n x CN (bool)
        # 1/q if q alleles has equal max probility
        belong_norm        = (belong / belong.sum(axis=2)[:, :, None])         # size: reads x top_n x CN (float)
        # sum across read and normalize by n
        allele_n_top_frac  = belong_norm.sum(axis=0) / self.log_probs.shape[0] # size: top_n x CN

        # sort and save the result
        self.result.append(TypingResult(
            n              = len(self.result) + 1,
            value          = allele_n_top_value,
            value_sum_indv = allele_n_top_sum,
            allele_id      = allele_n_top_id,
            allele_name    = self.mapAlleleIDs(allele_n_top_id),
            allele_prob    = allele_n_top_prob,
            fraction       = allele_n_top_frac,
        ).sortByScoreAndEveness())
        # breakpoint()
        return self.result[-1]

    def plot(self, title: str = "") -> list[go.Figure]:
        """ Plot the probility of read belonging to allele """
        probs = self.probs
        norm_probs = probs / probs.sum(axis=1, keepdims=True)
        # log_probs = np.log10(norm_probs)
        # plot prob
        figs = [
            px.imshow(np.log(probs), text_auto=True, aspect="auto", title=title),
            # px.imshow(norm_probs ** 0.1, text_auto=True, aspect="auto"),
            # px.imshow(np.log(norm_probs), text_auto=True, aspect="auto"),
        ]
        # figs.extend([dashbio.Clustergram(
        #     data=probs,
        #     column_labels=allele_names,
        #     hidden_labels='row',
        #     width=2000,
        #     height=1000,
        #     cluster="row",
        # )])
        return figs


class AlleleTypingExonFirst(AlleleTyping):
    only_exon_variant = False

    @staticmethod
    def aggrVariantsByAllele(variants: list[Variant]) -> dict[tuple, list[str]]:
        """ variant's alleles to dict[allele's variants, allele] """
        allele_variants = defaultdict(list)
        for variant in variants:
            for allele in variant.allele:
                allele_variants[allele].append(variant.id)
        variantset_to_allele = defaultdict(list)
        for allele, allele_vars in allele_variants.items():
            variantset_to_allele[tuple(sorted(set(allele_vars)))].append(allele)  # type: ignore
        return variantset_to_allele

    @staticmethod
    def removeIntronVariant(reads: list[PairRead], exon_variants: list[Variant]) -> list[PairRead]:
        """ Only exon variants will be preserved """
        id_map = {v.id: v for v in exon_variants}
        new_reads = copy.deepcopy(reads)
        for read in new_reads:
            read.lpv = [v for v in read.lpv if v in id_map]
            read.lnv = [v for v in read.lnv if v in id_map]
            read.rpv = [v for v in read.rpv if v in id_map]
            read.rnv = [v for v in read.rnv if v in id_map]
        return new_reads

    @staticmethod
    def createInverseMapping(allele_group: dict[str, list[str]]) -> dict[str, str]:
        """ {a: [b,c,d]} -> {b:a, c:a, d:a} """
        allele_map_to_group = {}
        for group, alleles in allele_group.items():
            for allele in alleles:
                allele_map_to_group[allele] = group
        return allele_map_to_group

    @staticmethod
    def removeDuplicateAllele(variants: list[Variant], allele_map: dict[str, str]) -> list[Variant]:
        """
        Some alleles are in the same group ({b:a, c:a}),
        so remove it in the variant.allele
        """
        variants = copy.deepcopy(variants)
        for variant in variants:
            variant.allele = list(set(filter(None, [allele_map.get(v, "") for v in variant.allele])))
        return variants

    def __init__(self, reads: list[PairRead], variants: list[Variant],
                 exon_only: bool = False,
                 candidate_set: str = "first_score",
                 candidate_set_threshold: float = 1.1):
        """ Extracting exon alleles """
        # extract exon variants
        exon_variants = [v for v in variants if v.in_exon]
        intron_variants = [v for v in variants if not v.in_exon]

        # cleanup reads (only exon variant preserved)
        exon_reads = self.removeIntronVariant(reads, exon_variants)
        exon_reads = self.removeEmptyReads(exon_reads)

        # aggr same alleles that has the same variantset
        variantset_to_allele = self.aggrVariantsByAllele(exon_variants)
        other_allele = self.collectAlleleNames(variants) \
                       - self.collectAlleleNames(exon_variants)
        if other_allele:  # special case
            variantset_to_allele[tuple()] = sorted(other_allele)
        self.allele_group = {alleles[0]: alleles for alleles in variantset_to_allele.values()}
        exon_variants = self.removeDuplicateAllele(variants, self.createInverseMapping(self.allele_group))
        from pprint import pprint
        pprint(self.allele_group)

        # same as before
        super().__init__(exon_reads, exon_variants)
        self.first_set_only = candidate_set == "first_score"
        if candidate_set == "same_score":
            self.candidate_set_max_score_ratio = 1.
        elif candidate_set == "max_score_ratio":
            self.candidate_set_max_score_ratio = candidate_set_threshold
        else:
            self.candidate_set_max_score_ratio = 1.1

        if not exon_only:
            self.full_model: AlleleTyping | None = AlleleTyping(reads, variants)
            self.full_model.top_n = 30  # TODO
        else:
            self.full_model = None

    def typingIntron(self, exon_candidates: list[list[str]], verbose=True) -> AlleleTyping:
        assert self.full_model
        model = copy.deepcopy(self.full_model)
        for cand in exon_candidates:
            res = model.addCandidate(cand)
            if verbose:
                res.print()
        return model

    def typing(self, cn: int) -> TypingResult:
        """
        Typing cn alleles.

        Returns:
          The top-n allele-set are best fit the reads (with maximum probility).
            Each set has CN alleles.
        """
        print("Exon:")
        result = super().typing(cn)
        result.fillNameGroup(self.allele_group)
        result.print()

        # run full typing as before but using selected alleles
        if self.full_model is None:
            # TODO: how to deal with exon-only sequences and
            # how to output (change selectBest() ? )
            return result
        assert cn == result.n
        if self.first_set_only:
            full_model = self.typingIntron(result.allele_name_group[0])
            self.result.extend(full_model.result)
            return self.result[-1]

        # all same-score KIR exon allele set
        max_score = result.value.max()
        mixed_result = []
        for i in range(len(result.value)):
            print(f"Typing Intron of set {i}")
            if result.value[i] < max_score * self.candidate_set_max_score_ratio:
                continue
            full_model = self.typingIntron(result.allele_name_group[i], verbose=False)
            self.result.extend(full_model.result)
            mixed_result.append(full_model.result[-1])
        print(f"{len(mixed_result)=}")

        # Merge result and sort it (same sorting method as addCandidate)
        result = TypingResult(
            n              = mixed_result[0].n,
            value          = np.concatenate([res.value          for res in mixed_result]),
            value_sum_indv = np.concatenate([res.value_sum_indv for res in mixed_result]),
            allele_id      = np.concatenate([res.allele_id      for res in mixed_result]),
            allele_name    = list(chain.from_iterable(res.allele_name for res in mixed_result)),
            allele_prob    = np.concatenate([res.allele_prob    for res in mixed_result], axis=1),
            fraction       = np.concatenate([res.fraction       for res in mixed_result]),
        )
        return result.sortByScoreAndEveness()
