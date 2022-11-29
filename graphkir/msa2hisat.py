"""
MSA -> hisat index
"""
import itertools
from typing import ClassVar, TextIO, Any
from dataclasses import dataclass, field

from Bio import SeqIO
from pyhlamsa import Genemsa
from .utils import runDocker
from .kir_msa import readFromMSAs


@dataclass
class Variant:
    """ Save variant's information """
    # basic
    pos: int  # position
    typ: str  # type of variant: insertion, deletion, or single(SNP)
    ref: str  # reference (KIR gene)
    val: None | int | str = None  # val: alt or deletion's length
    id: None | str = None  # variant ID
    length: int = 0  # for gK_hisat2

    # star-allele level
    allele: list[str] = field(default_factory=list)  # the star alleles has this variant
    freq: None | float = None  # allele frequency
    ignore: None | bool = False  # Ture if freq < min_freq_threshold
    min_freq_threshold: ClassVar[float] = 0.1
    in_exon: bool = False  # is the variant in exon

    # graph level
    count: ClassVar[int] = 0  # variant id
    haplo_id: ClassVar[int] = 0  # haplotype id
    novel_id: ClassVar[int] = 0  # novel allele id (increase when find new variant in in hisat2.py)

    # const
    order_type: ClassVar[dict[str, int]] = {"insertion": 0, "single": 1, "deletion": 2, "match": 3}
    order_nuc: ClassVar[dict[str, int]] = {"A": 0, "C": 1, "G": 2, "T": 3}

    def __lt__(self, othr: object) -> bool:
        if not isinstance(othr, Variant):
            return NotImplemented
        return (self.ref, self.pos, self.order_type[self.typ], self.val) \
             < (othr.ref, othr.pos, self.order_type[othr.typ], othr.val)

    def __eq__(self, othr: object) -> bool:
        if not isinstance(othr, Variant):
            return NotImplemented
        return (self.pos, self.ref, self.typ, self.val) \
            == (self.pos, self.ref, self.typ, self.val)

    def __hash__(self) -> int:
        return hash((self.pos, self.ref, self.typ, self.val))


def count2Freq(base_counts: list[list[int]]) -> list[dict[str, float]]:
    """
    Calculate per base frequency by its counts

    e.g. `[ [1,2,3,4, 0] ]` => `[ {'A': 0.1, 'T': 0.2, ...} ] `
    """
    freq = []
    for bc in base_counts:
        # normalize
        bc_norm = [c / sum(bc) for c in bc]
        # to dict (remove freq=0)
        freq.append(dict(filter(lambda i: i[1], zip("ATCG-", bc_norm))))
    return freq


def getVariantsFromSeqs(ref_seq: str, allele_seq: str, ref_name: str) -> list[Variant]:
    """ Call the variants of allele_seq from ref_seq """
    variants = []
    pre_v = None
    ref_skip_base = 0
    for i in range(len(ref_seq)):
        a = ref_seq[i]
        b = allele_seq[i]
        now_v = None
        if a != '-' and b != '-' and a != b:
            now_v = Variant(typ="single",    pos=i - ref_skip_base, val=b, ref=ref_name)
        elif a != '-' and b == '-':
            now_v = Variant(typ="deletion",  pos=i - ref_skip_base, val=1, ref=ref_name)
        elif a == '-' and b != '-':
            now_v = Variant(typ="insertion", pos=i - ref_skip_base, val=b, ref=ref_name)
        # a == '-' and b == '-'

        if a == '.':
            ref_skip_base += 1

        # merge pre_v and now_v
        if pre_v and now_v and pre_v.typ == now_v.typ and pre_v.typ != "single":
            pre_v.val += now_v.val
        # cannont merge -> previous variant done
        elif pre_v:
            variants.append(pre_v)
            pre_v = now_v
        else:
            pre_v = now_v
    if pre_v:  # this varinat is done
        variants.append(pre_v)

    return variants


def msa2Variant(msa: Genemsa) -> tuple[list[Variant], dict[str, list[Variant]]]:
    """
    MSA to variants

    This will return a tuple contains

    1. The list of variants (unique)
    2. A dictionary that key is allele name and value is its list of variants
    """
    ref_name, ref_seq = msa.get_reference()

    # get variantion from each allele
    variants_per_allele = {}
    for allele_name in msa.get_sequence_names():
        if allele_name == ref_name:
            continue
        variants_per_allele[allele_name] = \
            getVariantsFromSeqs(ref_seq, msa.get(allele_name), ref_name)

    # sort variantion and change it to dict
    # and add reference
    variants = sorted(set(itertools.chain.from_iterable(variants_per_allele.values())))
    variants_dict = {v: v for v in variants}

    # collect allele for the variant
    for allele, allele_variants in variants_per_allele.items():
        for v in allele_variants:
            variants_dict[v].allele.append(allele)

    # get frequency
    variants_dict = addFreqInVariant(msa, variants_dict)

    # add id
    for v in variants:
        variants_dict[v].id = f"hv{v.count}"
        Variant.count += 1

    # assign the same pointer for variant object
    for allele_name in variants_per_allele:
        variants_per_allele[allele_name] = \
            [variants_dict[v] for v in variants_per_allele[allele_name]]

    return list(variants_dict.values()), variants_per_allele


def addFreqInVariant(msa: Genemsa,
                     variants_dict: dict[Variant, Variant]
                     ) -> dict[Variant, Variant]:
    """
    Add allele frequency in variant (variants_dict).

    The frequency is calculated from MSA (msa).

    Note: Only SNPs are applied the minimal threshold
    """
    # calculate frequency to all bases
    base_counts = msa.calculate_frequency()
    base_freq = count2Freq(base_counts)

    # Get frequency for each variant
    for v in variants_dict:
        if v.typ == "deletion":
            base = '-'
        else:
            assert isinstance(v.val, (str))
            base = v.val[0]
        freq = base_freq[v.pos][base]
        variants_dict[v].freq = freq
        variants_dict[v].ignore = v.typ == "single" and freq < Variant.min_freq_threshold
    return variants_dict


def writeTSV(f: TextIO, data: list[Any]) -> None:
    """ Save list of data separated by tab """
    f.write('\t'.join(map(str, data)) + "\n")


def writeMsa(index_prefix: str, msa: Genemsa) -> None:
    """
    Write MSA-only information

    * `{index_prefix}_backbone.fa`
    * `{index_prefix}_sequences.fa`
    * `{index_prefix}.allele`
    * `{index_prefix}.partial`
    * `{index_prefix}.locus`
    """
    ref_name = msa.get_reference()[0]

    # Backbone Sequence
    with open(index_prefix + "_backbone.fa", 'a') as seq_backbone_f:
        SeqIO.write(msa.select_allele([ref_name]).to_records(gap=False),
                    seq_backbone_f, "fasta")

    # All sequence
    with open(index_prefix + "_sequences.fa", 'a') as seq_f:
        SeqIO.write(msa.copy().remove(ref_name).to_records(gap=False),
                    seq_f, "fasta")

    # Allele allele name
    with open(index_prefix + ".allele", 'a') as allele_name_f:
        allele_name_f.writelines(
            map(lambda i: i + "\n",  # type: ignore
                filter(lambda i: i != ref_name, msa.get_sequence_names()))
        )

    # Exon-only allele name
    with open(index_prefix + ".partial", 'a') as _:
        pass

    # Every msa has it's own reference
    with open(index_prefix + ".locus", 'a') as locus_f:
        # TODO: bounardy
        exon_pos = []
        for b in msa.blocks:
            pos = msa.get_block_position(b)
            if b.type == "exon":
                exon_pos.append((pos + 1, pos + 1 + b.length))
        exon_str = " ".join(f"{s}-{e}" for s, e in exon_pos)
        leng = msa.get_length()
        writeTSV(locus_f, [ref_name, ref_name, 0, leng, leng, exon_str, "+"])


def writeVariant(index_prefix: str, variants: list[Variant]) -> None:
    """
    Write variant information

    * `{index_prefix}.snp`
    * `{index_prefix}.snp.freq`
    * `{index_prefix}.link`
    """
    with open(index_prefix + ".snp", 'a') as snp_f:
        with open(index_prefix + ".index.snp", 'a') as snp_index_f:
            for v in variants:
                variant_tsv = [v.id, v.typ, v.ref, v.pos, v.val]
                if not v.ignore:
                    writeTSV(snp_index_f, variant_tsv)
                writeTSV(snp_f, variant_tsv)

    with open(index_prefix + ".snp.freq", 'a') as snp_freq_f:
        for v in variants:
            writeTSV(snp_freq_f, [v.id, f"{v.freq:.2f}"])

    with open(index_prefix + ".link", 'a') as snp_link_f:
        for v in variants:
            writeTSV(snp_link_f, [v.id, ' '.join(v.allele)])


def writeHaplo(index_prefix: str,
               variants_per_allele: dict[str, list[Variant]],
               msa: Genemsa) -> None:
    """
    Write star allele information

    * `{index_prefix}.haplotype`
    """
    with open(index_prefix + ".haplotype", 'a') as haplo_f:
        for variants in variants_per_allele.values():
            variants_index = [v for v in variants if not v.ignore]
            if variants_index:
                left = min([v.pos for v in variants_index])
                right = max([v.pos if v.typ != "deletion"
                             else v.pos + v.val - 1  # type: ignore
                             for v in variants_index])
                ref = variants[0].ref
            else:
                continue
                # skip if no variants in the allele
                # TODO: check this implemetation
                # old code: if no variants -> use 0 - msa_length
                # left = 0
                # right = msa.get_length() - 1
                # ref = msa.get_reference()[0]
            ids = ','.join(str(v.id) for v in variants_index)
            writeTSV(haplo_f, [f"ht{Variant.haplo_id}",
                               ref, left, right, ids])
            Variant.haplo_id += 1


def clearBeforeWrite(index_prefix: str) -> None:
    """ Clear the output file """
    extension = [".snp", ".index.snp", ".snp.freq", ".link", "_backbone.fa",
                 "_sequences.fa", ".allele", ".partial", ".locus", ".haplotype"]
    for ext in extension:
        with open(index_prefix + ext, 'w') as _:
            pass


def msa2HisatReference(msa_prefix: str, index_prefix: str) -> None:
    """
    Transfer MSA to hisat2 format

    Args:
      msa_prefix: The prefix of MSAs
      index_prefix: The prefix of output formats
    """
    Variant.min_freq_threshold = 0.1
    genes = readFromMSAs(msa_prefix)
    clearBeforeWrite(index_prefix)

    for gene, msa in genes.items():
        print("Reading", gene)

        # check
        for seq in msa.alleles.values():
            assert all(i in "ATCG-" for i in seq)
            assert len(seq) == msa.get_length()
        # Get reference and check
        ref_name, ref_seq = msa.get_reference()
        assert "BACKBONE" in ref_name
        assert all(i in "ATCG" for i in ref_seq)

        # main
        variants, variants_per_allele = msa2Variant(msa)

        # write to hisat format for hisat build
        writeMsa    (index_prefix, msa)
        writeVariant(index_prefix, variants)
        writeHaplo  (index_prefix, variants_per_allele, msa=msa)


def buildHisatIndex(name: str, output_name: str, threads: int = 1) -> None:
    """ Run hisat2-build, input and output are prefix of filenames """
    runDocker("hisat", f"""\
              hisat2-build {name}_backbone.fa \
                           --snp {name}.index.snp \
                           --haplotype {name}.haplotype \
                           -p {threads} --verbose \
                           {output_name} """)
