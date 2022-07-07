# The idea is copy from exon_vars in hisat-genotype/hisatgenotype_typing_process.py
# linnil1 2022/4/9
# I remove some features
# * Exon-related code
# * Find the position of target sequence in reference (and shift our position)
# * collapsed allele
# * reverse
# * split haplotype (break haplotype if two variant are too far away
import os
from glob import glob
import itertools
from typing import Union, ClassVar
from dataclasses import dataclass, field

from Bio import SeqIO
from pyhlamsa import msaio
from kg_utils import threads, runDocker, getSamples


@dataclass
class Variant:
    pos: int
    typ: str
    ref: str = None
    val: Union[int, str] = None

    allele: list = field(default_factory=list)
    id: str = None
    freq: float = None
    ignore: bool = False  # Ture if freq < min_freq_threshold

    order_type: ClassVar[dict] = {"insertion": 0, "single": 1, "deletion": 2}
    order_nuc: ClassVar[dict] = {"A": 0, "C": 1, "G": 2, "T": 3}
    min_freq_threshold: ClassVar[float] = 0.1
    count: ClassVar[int] = 0

    def __lt__(self, othr):
        return (self.ref, self.pos, self.order_type[self.typ], self.val) \
             < (othr.ref, othr.pos, self.order_type[othr.typ], othr.val)

    def __eq__(self, othr):
        return (self.pos, self.ref, self.typ, self.val) \
            == (self.pos, self.ref, self.typ, self.val)

    def __hash__(self):
        return hash((self.pos, self.ref, self.typ, self.val))


def count2Freq(base_counts):
    freq = []
    for bc in base_counts:
        # normalize
        bc = [c / sum(bc) for c in bc]
        # to dict (remove freq=0)
        freq.append(dict(filter(lambda i: i[1], zip("ATCG-", bc))))
    return freq


def getVariantsFromSeqs(ref_seq, allele_seq):
    """
    Compare two sequence and return a list of variants

    Variant canbe  'single', 'deletion', 'insertion'
    """
    variants = []
    pre_v = None
    ref_skip_base = 0
    for i in range(len(ref_seq)):
        a = ref_seq[i]
        b = allele_seq[i]
        now_v = None
        if a != '-' and b != '-' and a != b:
            now_v = Variant(typ="single", pos=i - ref_skip_base, val=b)
        elif a != '-' and b == '-':
            now_v = Variant(typ="deletion", pos=i - ref_skip_base, val=1)
        elif a == '-' and b != '-':
            now_v = Variant(typ="insertion", pos=i - ref_skip_base, val=b)
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


def msa2Variant(msa):
    """
    Main: Transfer msa to hisat2 format
    """
    ref_name, ref_seq = msa.get_reference()

    # get variantion from each allele
    variants_per_allele = {}
    for allele_name in msa.get_sequence_names():
        if allele_name == ref_name:
            continue
        variants_per_allele[allele_name] = \
            getVariantsFromSeqs(ref_seq, msa.get(allele_name))
        # add reference
        for v in variants_per_allele[allele_name]:
            v.ref = ref_name

    # sort variantion and change it to dict
    # and add reference
    variants = sorted(set(itertools.chain(*variants_per_allele.values())))
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

    return variants_dict, variants_per_allele


def addFreqInVariant(msa, variants_dict):
    """ Add alllele frequency in variant """
    # calculate frequency to all bases
    base_counts = msa.calculate_frequency()
    base_freq = count2Freq(base_counts)

    # Get frequency for each variant
    for v in variants_dict:
        if v.typ == "deletion":
            base = '-'
        else:
            base = v.val[0]
        freq = base_freq[v.pos][base]
        variants_dict[v].freq = freq
        variants_dict[v].ignore = v.typ == "single" and freq < Variant.min_freq_threshold
    return variants_dict


def writeTSV(f, data):
    f.write('\t'.join(map(str, data)) + "\n")


def writeMsa(index_out, msa):
    ref_name = msa.get_reference()[0]

    # Backbone Sequence
    with open(index_out + "_backbone.fa", 'a') as seq_backbone_f:
        SeqIO.write(msa.select_allele([ref_name]).to_records(gap=False),
                    seq_backbone_f, "fasta")

    # All sequence
    with open(index_out + "_sequences.fa", 'a') as seq_f:
        SeqIO.write(msa.copy().remove(ref_name).to_records(gap=False),
                    seq_f, "fasta")

    # Allele allele name
    with open(index_out + ".allele", 'a') as allele_name_f:
        allele_name_f.writelines(
            map(lambda i: i + "\n",
                filter(lambda i: i != ref_name, msa.get_sequence_names()))
        )

    # Exon-only allele name
    with open(index_out + ".partial", 'a') as partial_allele_name_f:
        pass

    # Every msa has it's own reference
    with open(index_out + ".locus", 'a') as locus_f:
        # TODO: bounardy
        exon_pos = []
        for b in msa.blocks:
            pos = msa.get_block_position(b)
            if b.type == "exon":
                exon_pos.append((pos + 1, pos + 1 + b.length))
        exon_str = " ".join(f"{s}-{e}" for s, e in exon_pos)
        leng = msa.get_length()
        writeTSV(locus_f, [ref_name, ref_name, 0, leng, leng, exon_str, "+"])


def writeVariant(index_out, variants_dict):
    with open(index_out + ".snp", 'a') as snp_f:
        with open(index_out + ".index.snp", 'a') as snp_index_f:
            for v in variants_dict.values():
                variant_tsv = [v.id, v.typ, v.ref, v.pos, v.val]
                if not v.ignore:
                    writeTSV(snp_index_f, variant_tsv)
                writeTSV(snp_f, variant_tsv)

    with open(index_out + ".snp.freq", 'a') as snp_freq_f:
        for v in variants_dict.values():
            writeTSV(snp_freq_f, [v.id, f"{v.freq:.2f}"])

    with open(index_out + ".link", 'a') as snp_link_f:
        for v in variants_dict.values():
            writeTSV(snp_link_f, [v.id, ' '.join(v.allele)])


def writeHaplo(index_out, variants_per_allele, msa):
    haplo_id = 0  # set variable in function (static)
    if hasattr(writeHaplo, "haplo_id"):
        haplo_id = writeHaplo.haplo_id

    with open(index_out + ".haplotype", 'a') as haplo_f:
        for allele, variants in variants_per_allele.items():
            variants_index = [v for v in variants if not v.ignore]
            if variants_index:
                left = min([v.pos for v in variants_index])
                right = max([v.pos if v.typ != "deletion" else v.pos + v.val - 1 for v in variants_index])
                ref = variants[0].ref
            else:
                continue
                left = 0
                right = msa.get_length() - 1
                ref = msa.get_reference()[0]
            ids = ','.join([v.id for v in variants_index])
            writeTSV(haplo_f, [f"ht{haplo_id}",
                               ref, left, right, ids])
            haplo_id += 1
    writeHaplo.haplo_id = haplo_id


def buildHisat(name):
    runDocker("hisat", f"""\
              hisat2-build {name}_backbone.fa \
                           --snp {name}.index.snp \
                           --haplotype {name}.haplotype \
                           -p {threads} --verbose \
                           {name}.graph """)


def main(index, force=False):
    # remove `.save`
    # relate to min_freq_threshold
    index_out = index[:-5] + ".mut01"
    Variant.min_freq_threshold = 0.1

    # remove old file, because I write file by append not overwrite
    if not force:
        if getSamples(index_out + ".graph", strict=False):
            return index_out

    for i in getSamples(index_out, strict=False):
        os.remove(i)
    for i in glob(index_out + "_*"):
        os.remove(i)

    # get all genes
    # genes = sorted([f.split('.')[-3] for f in glob(index + ".*.save.json")])
    genes = sorted([i[1] for i in getSamples(index, ".json", return_name=True)])
    print(genes)

    for gene in genes:
        # read
        print("Reading", gene)
        msa = msaio.load_msa(f"{index}.{gene}.fa", f"{index}.{gene}.json")

        # check
        for seq in msa.alleles.values():
            assert all(i in "ATCG-" for i in seq)
            assert len(seq) == msa.get_length()

        # Get reference and check
        ref_name, ref_seq = msa.get_reference()
        assert "BACKBONE" in ref_name
        assert all(i in "ATCG" for i in ref_seq)

        # main
        variants_dict, variants_per_allele = msa2Variant(msa)

        # write to hisat format for hisat build
        writeMsa(index_out, msa)
        writeVariant(index_out, variants_dict)
        writeHaplo(index_out, variants_per_allele, msa=msa)

    # hisat build
    buildHisat(index_out)
    return index_out


if __name__ == "__main__":
    index = "index1/kir_2100_2dl1s1.save"
    index = main(index)
    print(index)
