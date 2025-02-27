"""
Typing the sample's alleles using variants file and CN file as input
"""
from typing import Any
from collections import defaultdict
import json

from .utils import NumpyEncoder, logger
from .msa2hisat import Variant
from .hisat2 import loadReadsAndVariantsData, removeMultipleMapped, PairRead
from .typing_mulit_allele import AlleleTyping, AlleleTypingExonFirst
from .typing_em import preprocessHisatReads, hisat2TypingPerGene, printHisatTyping


def groupReads(reads: list[PairRead]) -> dict[str, list[PairRead]]:
    """Group the reads by reference name (i.e. group by gene)"""
    gene_reads = defaultdict(list)
    for read in reads:
        gene_reads[read.backbone].append(read)
    return gene_reads


def groupVariants(variants: list[Variant]) -> dict[str, list[Variant]]:
    """Group the reads by reference name (i.e. group by gene)"""
    gene_variants = defaultdict(list)
    for variant in variants:
        gene_variants[variant.ref].append(variant)
    return gene_variants


class Typing:
    """Abstract class for typing allele"""

    def __init__(self) -> None:
        """Read sample's variants or allele abundance result"""
        self._result: dict[str, Any] = {}

    def typingPerGene(self, gene: str, cn: int) -> tuple[list[str], int]:
        """Typing for allele within gene with given CN"""
        raise NotImplementedError

    def typing(self, gene_cn: dict[str, int], min_reads_num: int = 100) -> tuple[list[str], list[str]]:
        """
        Allele Typing for all genes by given CN value per gene (`gene_cn`).

        Returns:
          list of alleles in this sample
        """
        predict_alleles = []
        warning_genes = []
        for gene, cn in gene_cn.items():
            if not cn:
                continue
            # debug
            # if "KIR3DP1" not in gene:
            #     continue
            predict_allele, reads_num = self.typingPerGene(gene, cn)
            predict_alleles.extend(predict_allele)
            if reads_num < min_reads_num:  # 100(reads) * 300(base/read) ~= 5x (min 5k length)
                warning_genes.append(gene)

        return predict_alleles, warning_genes

    def save(self, filename: str) -> None:
        """Save data in file"""
        with open(filename, "w") as f:
            json.dump(self._result, f, cls=NumpyEncoder)

    def getAllPossibleTyping(self) -> list[dict[Any, Any]]:
        """
        Return all possible set of allele typing.
        Call typing() before calling this function
        """
        raise NotImplementedError


class TypingWithPosNegAllele(Typing):
    """Our proposed allele typing method"""

    def __init__(
        self,
        filename_variant_json: str,
        top_n: int = 300,
        multiple: bool = False,
        exon_first: bool = False,
        exon_only: bool = False,
        exon_candidate_threshold: float = .9,
        variant_correction: bool = False,
    ):
        """Read all reads and variants from the json file (From .hisat2.py)"""
        super().__init__()
        reads_data = loadReadsAndVariantsData(filename_variant_json)
        if not multiple:
            reads_data = removeMultipleMapped(reads_data)
        self._top_n = top_n
        self._gene_reads = groupReads(reads_data["reads"])
        self._gene_variants = groupVariants(reads_data["variants"])
        self._exon_first = exon_first
        self._exon_only = exon_only
        self._exon_candidate_threshold = exon_candidate_threshold
        self._variant_correction = variant_correction

    def typingPerGene(self, gene: str, cn: int) -> tuple[list[str], int]:
        """Select reads belonged to the gene and typing it"""
        logger.debug(f"[Allele] {gene=} {cn=}")
        if not self._exon_first and not self._exon_only:
            typ = AlleleTyping(
                self._gene_reads[gene],
                self._gene_variants[gene],
                gene,
                cn,
                top_n=self._top_n,
                variant_correction=self._variant_correction,
            )
            if not typ.homo:
                res = typ.typing(cn)
                self._result[gene] = typ.result
                alleles = res.selectBest()
                # KIR2DL1*BACKBONE -> KIR2DL1
                alleles = [i if i != "fail" else f"{pure_gene}*" for i in alleles]
            else:
                alleles = [f"{typ.typingHomo()}" for i in range(cn)]
        else:
            typ = AlleleTypingExonFirst(
                self._gene_reads[gene],
                self._gene_variants[gene],
                top_n=self._top_n,
                exon_only=self._exon_only,
                candidate_set_threshold=self._exon_candidate_threshold,
            )
            res = typ.typing(cn)
            self._result[gene] = typ.result
            # return res.selectBest(filter_minor=True)
            alleles = res.selectBest()
            # KIR2DL1*BACKBONE -> KIR2DL1
            pure_gene = gene.split("*")[0]
            alleles = [i if i != "fail" else f"{pure_gene}*" for i in alleles]
        return alleles, typ.getReadsNum()

    def getAllPossibleTyping(self) -> list[dict[Any, Any]]:
        """
        Return all possible set of allele typing.
        Call typing() before calling this function
        """
        possible_list = []
        for gene, result in self._result.items():
            for rank, (value, alleles) in enumerate(result[-1].selectAllPossible(.9)):
                simple_result = {
                    "gene": gene,
                    "rank":  rank,
                    "value": value,
                }
                for i, allele in enumerate(alleles):
                    simple_result[str(i + 1)] = allele
                possible_list.append(simple_result)
        return possible_list


class TypingWithReport(Typing):
    """Typing alleles from report by abundance value"""

    def __init__(self, filename_variant_json: str):
        """Read report files (json)"""
        super().__init__()
        reads_data = loadReadsAndVariantsData(filename_variant_json)
        reads_data = removeMultipleMapped(reads_data)
        self._gene_reads = preprocessHisatReads(reads_data)

    def typingPerGene(self, gene: str, cn: int) -> tuple[list[str], int]:
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
        return called_alleles, len(self._gene_reads[gene])

    def save(self, filename: str) -> None:
        """save additional report txt"""
        super().save(filename)
        name = filename
        if filename.endswith(".json"):
            name = filename[:-5]
        with open(name + ".txt", "w") as f:
            printHisatTyping(self._result, file=f)


def selectKirTypingModel(
    method: str,
    filename_variant_json: str,
    **kwargs: Any,
) -> Typing:
    """Select and Init typing model"""
    if method == "full":
        return TypingWithPosNegAllele(filename_variant_json, **kwargs)
    if method.startswith("exonfirst"):
        fields = method.split("_")
        threshold = 0.0
        if len(fields) == 2:  # e.g. exonfirst_1.2
            threshold = float(method[len("exonfirst_") :])
        return TypingWithPosNegAllele(
            filename_variant_json,
            exon_first=True,
            exon_candidate_threshold=threshold,
            **kwargs,
        )
    if method == "em":
        return TypingWithReport(filename_variant_json)
    raise NotImplementedError
