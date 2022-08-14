"""
Our main typing method
"""
from itertools import chain
from dataclasses import dataclass
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from gk_hisat2 import PairRead, Variant


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
    """
    n: int
    value:       np.ndarray       # size: top_n
    allele_id:   np.ndarray       # size: top_n x n
    allele_name: list[list[str]]  # size: top_n x n
    allele_prob: np.ndarray       # size: read x top_n
    fraction:    np.ndarray       # size: top_n x n

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
              id   1 name KIR3DS1*0130108      fraction 0.35294117647058826
              id   2 name KIR3DS1*078          fraction 0.6470588235294118
            1 -52.37587690948907
              id   0 name KIR3DS1*0130107      fraction 0.4215686274509804
              id   2 name KIR3DS1*078          fraction 0.5784313725490197
            ```
        """
        print(self.n)
        top_n = len(self.value)
        n = self.allele_id.shape[1]

        for rank in range(min(top_n, num)):
            print("Rank", rank, "probility", self.value[rank])

            for i in range(n):
                print(f"  id {self.allele_id[rank][i]:3}"
                      f"  name {self.allele_name[rank][i]:20s}"
                      f"  fraction {self.fraction[rank][i]:.5f}")
                # f" unique_count {unique_count}")


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

        allele_names = self.collectAlleleNames(reads)
        self.id_to_allele: dict[int, str] = dict(enumerate(sorted(allele_names)))
        self.allele_to_id: dict[str, int] = {j: i for i, j in self.id_to_allele.items()}

        self.probs = self.reads2AlleleProb(reads)
        self.log_probs = np.log10(self.probs)

        self.result: list[TypingResult] = []
        # , allele_length: dict[str, float]):
        # allele_length = np.array([allele_length_map[self.allele_name_map_rev[i]]
        # for i in range(len(self.allele_name_map_rev))])
        # self.allele_length = allele_length / 10000.  # just a non-important normalize

    def collectAlleleNames(self, reads: list[PairRead]) -> set[str]:
        """ Get all allele names from all the reads """
        allele_names = set()
        for read in reads:
            variants = chain(read.lpv, read.rpv, read.lnv, read.rnv)
            alleles = map(lambda i: self.variants[i].allele, variants)
            allele_names.update(set(chain.from_iterable(alleles)))
        return allele_names

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
    def argSortRow(data: np.ndarray) -> list[int]:
        """ Argsort the data by row """
        return sorted(range(len(data)), key=lambda i: tuple(data[i]))

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

    def addCandidate(self) -> TypingResult:
        """
        The step of finding the maximum likelihood allele-set that can fit all reads

        Once call this function, the number alleles in the allele set increate 1.
        """
        if not len(self.result):
            # first time: i.e. CN = 1
            allele_1_prob      = self.log_probs.sum(axis=0)                     # size: allele
            allele_1_top_index = np.argsort(allele_1_prob)[-self.top_n:][::-1]  # size: top_n
            allele_1_top_id    = np.array([[i] for i in allele_1_top_index])    # size: top_n x 1
            allele_1_top_prob  = self.log_probs[:, allele_1_top_index]          # size: reads x top_n

            self.result.append(TypingResult(
                n              = 1,
                value          = allele_1_prob[allele_1_top_index],             # size: top_n
                allele_id      = allele_1_top_id,
                allele_name    = self.mapAlleleIDs(allele_1_top_id),            # size: top_n x 1
                allele_prob    = allele_1_top_prob,
                fraction       = np.ones(allele_1_top_id.shape),                # size:  top_n x 1
            ))
            return self.result[-1]

        # get previous top-n allele
        allele_n_1_prob = self.result[-1].allele_prob  # size: read x top_n
        allele_n_1_id   = self.result[-1].allele_id    # size: top_n x CN

        # Maximum
        # (top-n x (allele-1 ... allele-n-1)) x (allele_1 ... allele_m)
        allele_n_prob = np.maximum(self.log_probs, allele_n_1_prob.T[:, :, None])  # size: top_n x read x allele
        allele_n_prob = allele_n_prob.sum(axis=1).flatten()                        # size: (top_n * allele)

        # Mean
        # allele_n = (log_probs + prob_1_top.T[:, :, None] * (prob_iter - 1)).sum(axis=1).flatten() / prob_iter

        # create id (size: top_n * allele x n) of the allele_n_prob
        n_allele = self.log_probs.shape[1]
        allele_n_id = np.hstack([
            # [[1],[3],[5]] -> [[1],[3],[5],[1],[3],[5]]
            np.repeat(allele_n_1_id, n_allele, axis=0),
            # [1,3,5] -> [1,3,5,1,3,5]
            np.tile(np.arange(n_allele), len(allele_n_1_id))[:, None],
        ])

        # unique the allele set
        unique_id_index    = self.uniqueAllele(allele_n_id)
        allele_n_id        = allele_n_id[unique_id_index]
        allele_n_prob      = allele_n_prob[unique_id_index]

        # top-n result
        allele_n_top_index = np.argsort(allele_n_prob)[-self.top_n:][::-1]   # size: top_n
        allele_n_top_id    = allele_n_id[allele_n_top_index]                 # size: top_n x n
        allele_n_top_prob  = self.log_probs[:, allele_n_top_id].max(axis=2)  # read x top_n x n -> size: read x top-n
        allele_n_top_value = allele_n_prob[allele_n_top_index]               # size: top_n

        # calculate fraction
        belong                = np.equal(self.log_probs[:, allele_n_top_id],
                                         allele_n_top_prob[:, :, None])            # size: reads x top-n x n (bool)
        # 1/q if q alleles has equal max probility
        belong_norm           = (belong / belong.sum(axis=2)[:, :, None])          # size: reads x top_n x n (float)
        # sum across read and normalize by n
        allele_n_top_fraction = belong_norm.sum(axis=0) / self.log_probs.shape[0]  # size: top_n x n

        # sorted by likelihood and the difference value toward evenly-distributed
        fraction_diff = allele_n_top_fraction - allele_n_top_fraction.mean(axis=1, keepdims=True)  # size: top_n x n
        fraction_diff = np.abs(fraction_diff).sum(axis=1)                                          # size: top_n
        rank_index = self.argSortRow(np.array([-allele_n_top_value, fraction_diff]).T)             # size: top_n

        # save result
        self.result.append(TypingResult(
            n=len(self.result) + 1,
            value=allele_n_top_value[rank_index],
            allele_id=allele_n_top_id[rank_index],
            allele_name=self.mapAlleleIDs(allele_n_top_id[rank_index]),
            allele_prob=allele_n_top_prob[:, rank_index],
            fraction=allele_n_top_fraction[rank_index],
        ))
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
