"""
Typing the sample's alleles using variants file and CN file as input
"""
import json
from typing import Any
from collections import defaultdict
from dataclasses import asdict

from .gk_utils import NumpyEncoder
from .gk_hisat2 import loadReadsAndVariantsData, PairRead, removeMultipleMapped
from .gk_multi_allele_typing import AlleleTyping
from .gk_hisat2em import preprocessHisatReads, hisat2TypingPerGene, printHisatTyping


def groupReads(reads: list[PairRead]) -> dict[str, list[PairRead]]:
    """ Group the reads by reference name (i.e. group by gene) """
    gene_reads = defaultdict(list)
    for read in reads:
        gene_reads[read.backbone].append(read)
    return gene_reads


class Typing:
    """ Abstract class for typing allele """
    def __init__(self):
        """ Read sample's variants or allele abundance result """
        self._result: dict[str, Any] = {}

    def typingPerGene(self, gene: str, cn: int) -> list[str]:
        """ Typing for allele within gene with given CN """
        raise NotImplementedError

    def typing(self, gene_cn: dict[str, int]) -> list[str]:
        """
        Allele Typing for all genes by given CN value per gene (`gene_cn`).

        Returns:
          list of alleles in this sample
        """
        predict_alleles = []
        for gene, cn in gene_cn.items():
            if not cn:
                continue
            predict_alleles.extend(self.typingPerGene(gene, cn))
        return predict_alleles

    def save(self, filename: str):
        """ Save data in file """
        with open(filename, "w") as f:
            json.dump(self._result, f, cls=NumpyEncoder)


class TypingWithPosNegAllele(Typing):
    """ Our proposed allele typing method """

    def __init__(self, filename_variant_json: str, multiple=False):
        """ Read all reads and variants from the json file (From gk_hisat2.py) """
        super().__init__()
        reads_data = loadReadsAndVariantsData(filename_variant_json)
        if not multiple:
            reads_data = removeMultipleMapped(reads_data)
        self._gene_reads = groupReads(reads_data['reads'])
        self._variants = reads_data['variants']

    def typingPerGene(self, gene: str, cn: int) -> list[str]:
        """ Select reads belonged to the gene and typing it """
        print(gene, cn)
        typ = AlleleTyping(self._gene_reads[gene], self._variants)
        res = typ.typing(cn)
        self._result[gene] = typ.result
        return res.selectBest()


class TypingWithReport(Typing):
    """ Typing alleles from report by abundance value """

    def __init__(self, filename_variant_json: str):
        """ Read report files (json) """
        super().__init__()
        # break down gk_hisat2_em
        reads_data = loadReadsAndVariantsData(filename_variant_json)
        reads_data = removeMultipleMapped(reads_data)
        self._gene_reads = preprocessHisatReads(reads_data)

    def typingPerGene(self, gene: str, cn: int) -> list[str]:
        """
        Typing for allele within gene with given CN

        Example:
          * case1: a1=0.66 a2=0.33
              * CN=1: a1
              * CN=2: a1 a1
              * CN=3: a1 a1 a2
          * case2: a1=0.90 a2=0.10
              * CN=1: a1
              * CN=2: a1 a1
              * CN=3: a1 a1 a1
        """
        # main typing
        report = hisat2TypingPerGene(self._gene_reads[gene])
        report = sorted(report, key=lambda i: -i.prob)

        # main calling logic wrote here
        est_prob = 1 / cn
        called_alleles = []
        for allele in report:
            pred_count = max(1, round(allele.prob / est_prob))
            for _ in range(min(cn, pred_count)):
                called_alleles.append(allele.allele)
            allele.cn = pred_count

            cn -= pred_count
            if cn <= 0:
                break

        self._result[gene] = report
        return called_alleles

    def save(self, filename: str):
        """ save additional report txt """
        super().save(filename)
        name = filename
        if filename.endswith(".json"):
            name = filename[:-5]
        with open(name + ".txt", "w") as f:
            printHisatTyping(self._result, file=f)


def selectKirTypingModel(method: str, filename_variant_json: str) -> Typing:
    """ Select and Init typing model """
    if method == "pv":
        return TypingWithPosNegAllele(filename_variant_json)
    elif method == "em":
        return TypingWithReport(filename_variant_json)
    raise NotImplementedError