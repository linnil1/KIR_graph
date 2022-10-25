from pprint import pprint
from typing import DefaultDict, Callable, Iterable
from itertools import chain
from functools import partial
from collections import defaultdict, Counter
from pysam import AlignmentFile
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyhlamsa import msaio, Genemsa
from graphkir.utils import samtobam
from graphkir.hisat2 import PairRead, Variant
from graphkir.kir_typing import TypingWithPosNegAllele
from graphkir.typing_mulit_allele import AlleleTyping
from kg_utils import runDocker
from kg_eval import readPredictResult


GroupPairRead = dict[tuple[str, ...], list[PairRead]]


def groupReadByAllele(
    typ: AlleleTyping, predict_alleles: list[str], reads: list[PairRead]
) -> GroupPairRead:
    allele_names = []
    allele_ids = []
    for i in predict_alleles:
        if i in typ.allele_to_id:
            allele_names.append(i)
            allele_ids.append(typ.allele_to_id[i])
    if not len(allele_names):
        return {}
    is_max = np.equal(
        typ.probs[:, allele_ids], typ.probs[:, allele_ids].max(axis=1)[:, None]
    )
    assign_reads = defaultdict(list)
    for read, max_id in zip(reads, is_max):
        selected_alleles = tuple(sorted(np.array(allele_names)[max_id]))
        assign_reads[selected_alleles].append(read)
    return assign_reads


def variantConfusionInRead(
    read: PairRead, allele: str, variants
) -> dict[str, list[str]]:
    variant_confusion: dict[str, list[str]] = {
        "novel": [],
        "tp": [],
        "tn": [],
        "fp": [],
        "fn": [],
    }
    for v in chain(read.lpv, read.rpv):
        if v.startswith("nv"):
            variant_confusion["novel"].append(v)
            continue
        if allele in variants[v].allele:
            variant_confusion["tp"].append(v)
        else:
            variant_confusion["fp"].append(v)
    for v in chain(read.lnv, read.rnv):
        if v.startswith("nv"):
            variant_confusion["novel"].append(v)
            continue
        if allele in variants[v].allele:
            variant_confusion["fn"].append(v)
        else:
            variant_confusion["tn"].append(v)
    return variant_confusion


def analyzeConfusion(
    allele: str, reads: list[PairRead], variants: dict[str, Variant]
) -> dict[str, int]:
    # stat of correntness count
    count = {
        "total": 0,
        "novel": 0,
        "tp": 0,
        "tn": 0,
        "fp": 0,
        "fn": 0,
    }
    for read in reads:
        variant_confusion = variantConfusionInRead(read, allele, variants)
        for stat, vs in variant_confusion.items():
            count[stat] += len(vs)
        count["total"] = count["tp"] + count["tn"] + count["fp"] + count["fn"]
    return count


def extractNovelVariants(
    allele: str, reads: list[PairRead], variants: dict[str, Variant], threshold: int = 4
) -> dict[str, dict[Variant, int]]:
    # fp/fn/novel allele
    variant_count: dict[str, DefaultDict[str, int]] = {
        "fp": defaultdict(int),
        "fn": defaultdict(int),
        "novel": defaultdict(int),
    }
    for read in reads:
        variant_confusion = variantConfusionInRead(read, allele, variants)
        for stat, vs in variant_confusion.items():
            if stat not in ["fp", "fn", "novel"]:
                continue
            for v in vs:
                variant_count[stat][v] += 1

    novel_variants = {}
    for stat, v_count in variant_count.items():
        selected_variants = {
            variants[v]: num for v, num in v_count.items() if num > threshold
        }
        novel_variants[stat] = selected_variants
    return novel_variants


def applyNovelVariant(
    allele: str,
    reads: list[PairRead],
    novel_variants: dict[str, dict[Variant, int]],
    msa: Genemsa,
    queryPileup: Callable,
) -> tuple[bool, str]:
    # sequences
    backbone_seq = msa.get_reference()[1]
    allele_seq_msa = msa.get(allele).replace("E", "")
    allele_seq = allele_seq_msa.replace("-", "")
    read_num = len(reads)
    print("  leng", len(allele_seq))
    print("  read", read_num)
    print("  avg depth", read_num * 500 / len(allele_seq))

    # probs = []
    # if "KIR3DP1" in read.l_sam:
    #     continue
    # if stat['fn'] or stat['fp']:
    #     print(read.l_sam.split("\t")[0])
    #     prob = typ.reads2AlleleProb([read])[0, [7, 33, 64]]
    #     print(prob)
    #     probs.append(prob)
    # print("sum", np.sum(np.log10(probs), axis=0))

    # apply variant
    is_novel = False
    for stat, v_count in novel_variants.items():
        print(" ", stat)
        if not v_count:
            print(
                "    []",
            )

        for v, num in v_count.items():
            print(f"    Apply {v.ref}:{v.pos} {v.typ}{v.val} id={v.id} count={num}")
            pos = v.pos
            print("      Reference", backbone_seq[pos - 5 : pos + 6])
            print("      Before:  ", allele_seq_msa[pos - 5 : pos + 6])

            # pileup
            pileup_read_base = queryPileup(pos=pos)
            selected_read_id = set(read.l_sam.split("\t", 1)[0] for read in reads)
            selected_read_base = {
                id: base
                for id, base in pileup_read_base.items()
                if id in selected_read_id
            }
            if not len(selected_read_base):
                print("      Ignore (depth = 0)")
                continue
            print(f"      Pileup depth {len(selected_read_base)}")
            base_count = Counter(selected_read_base.values())
            print(f"      Pileup base= {base_count}")
            max_count = max(base_count.values())
            max_bases = set(
                base for base, count in base_count.items() if count == max_count
            )

            if stat == "fp":
                base_ref = allele_seq_msa[pos]
                base_alt = v.val
                assert base_ref == allele_seq_msa[pos]
            elif stat == "fn":
                base_ref = allele_seq_msa[pos]
                base_alt = backbone_seq[pos]
                assert base_ref == allele_seq_msa[pos]
            elif stat == "novel":
                base_ref = allele_seq_msa[pos]
                base_alt = v.val
            if base_alt not in max_bases:
                print("      Ignore (base is not at highest abundance)")
                continue
            if v.typ != "single":
                print("      Skip")
                continue
            assert base_ref != base_alt
            assert base_alt
            assert type(base_alt) is str
            is_novel = True
            allele_seq_msa = allele_seq_msa[:pos] + base_alt + allele_seq_msa[pos + 1 :]
            print("      After:   ", allele_seq_msa[pos - 5 : pos + 6])
    print("  novel", is_novel)
    return is_novel, allele_seq_msa


def groupReadToBam(
    input_bam: str, output_bam: str, assign_reads: GroupPairRead
) -> None:
    output_name, bam_ext = output_bam.rsplit(".", 1)
    assert bam_ext == "bam"
    runDocker("samtools", f"samtools view -H {input_bam} -o {output_name}.sam")
    with open(f"{output_name}.sam", "a") as f:
        for alleles in assign_reads.keys():
            alleles_str = ",".join(alleles)
            f.write(f"@RG\tID:{alleles_str}\n")
        for alleles, reads in assign_reads.items():
            alleles_str = ",".join(alleles)
            for read in reads:
                f.write(read.l_sam + f"\tRG:Z:{alleles_str}\n")
                f.write(read.r_sam + f"\tRG:Z:{alleles_str}\n")
    samtobam(output_name)


def queryPileup(bamfile: AlignmentFile, ref: str, pos: int) -> dict[str, str]:
    base = {}
    for pileupcolumn in bamfile.pileup(ref, pos, pos + 1):
        if pileupcolumn.pos != pos:  # type: ignore
            continue
        for pileupread in pileupcolumn.pileups:  # type: ignore
            if not pileupread.is_del and not pileupread.is_refskip:
                base[
                    pileupread.alignment.query_name
                ] = pileupread.alignment.query_sequence[pileupread.query_position]
    # for i, j in base.items():
    #     print(f"{i:30s} {j}")
    return base


def writeFasta(sequences: Iterable[tuple[str, str]], output_fasta: str) -> None:
    SeqIO.write(
        [SeqRecord(Seq(seq), id=id, description="") for id, seq in sequences],
        output_fasta,
        "fasta",
    )


def typingNovel(
    input_name: str,
    msa_name: str,
    result_name: str,
    output_name: str,
) -> None:
    suffix_bam = ".no_multi"
    # read tsv
    predict_alleles = list(readPredictResult(result_name + ".tsv").values())[0]
    print(predict_alleles)
    # predict_alleles = [i for i in predict_alleles if i.startswith("KIR2DL5")]
    # predict_alleles = ["KIR2DL5A*00103", "KIR2DL5A*031"]

    # read json and bam
    bamfile = AlignmentFile(input_name + suffix_bam + ".bam", "rb")
    data = TypingWithPosNegAllele(input_name + ".json")

    assign_reads_all = {}
    allele_seqs = []
    for gene, reads in data._gene_reads.items():
        # msa
        msa_gene_name = msa_name + "." + gene.split("*")[0]
        msa = msaio.load_msa(msa_gene_name + ".fa", msa_gene_name + ".json")

        # read variants and group them by called allele
        typ = AlleleTyping(reads, data._gene_variants[gene], no_empty=False)
        assert typ.probs.shape[0] == len(reads)
        assign_reads = groupReadByAllele(typ, predict_alleles, reads)
        assign_reads_all.update(assign_reads)

        # reads per allele
        for alleles, reads in assign_reads.items():
            # the reads cannot be group into one of the allele
            if len(alleles) > 1:
                continue
            allele = alleles[0]
            print(allele)

            # staticstic
            print(" ", analyzeConfusion(allele, reads, typ.variants))
            # find fp/fn/novel variants
            novel_variants = extractNovelVariants(
                allele, reads, typ.variants, threshold=4
            )
            # fp/fn/novel variants -> seq
            is_novel, seq = applyNovelVariant(
                allele,
                reads,
                novel_variants,
                msa,
                partial(queryPileup, bamfile=bamfile, ref=gene),
            )
            # save seq
            allele_name = allele
            if is_novel:
                allele_name += "_new"
            allele_seqs.append((allele_name, seq))

    groupReadToBam(
        input_name + suffix_bam + ".bam", output_name + ".bam", assign_reads_all
    )
    writeFasta(allele_seqs, output_name + ".fa")


if __name__ == "__main__":
    ref_name = "index5/kir_2100_ab_2dl1s1.leftalign"
    input_name = "data5/linnil1_syn_30x_seed87.02.index5_kir_2100_ab_2dl1s1.leftalign.mut01.graph.variant.noerrcorr"
    suffix_called_alleles = ".no_multi.depth.p75.CNgroup_assume3DL3.pv.compare_sum"
    output_name = "tmp"
    typingNovel(input_name, ref_name, input_name + suffix_called_alleles, output_name)
