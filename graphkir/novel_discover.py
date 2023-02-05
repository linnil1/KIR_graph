"""
Discover the novel variant from read and called alleles

1. Assign reads to alleles
2. Find unmatched variants in those reads
3. print the stat
4. Try to apply the variant
"""
import sys
from typing import Iterable, TypedDict, TextIO
from itertools import chain
from collections import defaultdict, Counter

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pysam import AlignmentFile
import numpy as np
import pandas as pd

from pyhlamsa import Genemsa
from .utils import samtobam, runDocker
from .hisat2 import PairRead
from .msa2hisat import Variant
from .kir_typing import TypingWithPosNegAllele
from .typing_mulit_allele import AlleleTyping


GroupPairRead = dict[tuple[str, ...], list[PairRead]]


class NovelVariant(TypedDict):
    gene: str
    allele: str
    allele_count: int  # The i'th allele in the gene
    type: str
    variant: Variant
    pos: int
    count: int
    skip: bool
    skip_reason: str
    base_ref: str
    base_alt: str
    pileup: dict[str, int]


def groupReadByAllele(
    typ: AlleleTyping, predict_alleles: list[str], reads: list[PairRead]
) -> GroupPairRead:
    """
    Assign the reads by predict_alleles
    with probability calculated from allele typing
    """
    allele_names = []
    allele_ids = []
    for i in predict_alleles:
        if i in typ.allele_to_id:
            allele_names.append(i)
            allele_ids.append(typ.allele_to_id[i])
    if not allele_names:
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
    read: PairRead,
    allele: str,
    variants: dict[str, Variant],
) -> dict[str, list[str]]:
    """Find mismatch variants against assigned allele"""
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
            variant_confusion["tp"   ].append(v)
        else:
            variant_confusion["fp"   ].append(v)
    for v in chain(read.lnv, read.rnv):
        if v.startswith("nv"):
            variant_confusion["novel"].append(v)
            continue
        if allele in variants[v].allele:
            variant_confusion["fn"   ].append(v)
        else:
            variant_confusion["tn"   ].append(v)
    return variant_confusion


def statNovelConfusion(
    allele: str, reads: list[PairRead], variants: dict[str, Variant]
) -> dict[str, int]:
    """Stat mismatch variants against assigned allele"""
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


def extractNovelVariant(
    allele: str, reads: list[PairRead], variants: dict[str, Variant]
) -> dict[str, dict[Variant, int]]:
    """Get novel variant info from mismatch variants"""
    novel_variants: dict[str, list[Variant]] = {
        "novel": [],
        "fp":    [],
        "fn":    [],
    }
    for read in reads:
        variant_confusion = variantConfusionInRead(read, allele, variants)
        for stat, vs in variant_confusion.items():
            if stat not in novel_variants:
                continue
            novel_variants[stat] += [variants[v] for v in vs]

    count = {}
    for stat, variants_list in novel_variants.items():
        count[stat] = dict(Counter(variants_list))
    return count


def updateBaseRefAlt(
    variant: NovelVariant, backbone_seq: str, allele_seq: str
) -> NovelVariant:
    """Find ref/alt of the variant"""
    # determine ref and alt
    v = variant["variant"]
    base_ref = allele_seq[v.pos]
    if   variant["type"] == "fp":
        base_alt = v.val
    elif variant["type"] == "novel":
        base_alt = v.val
    elif variant["type"] == "fn":
        base_alt = backbone_seq[v.pos]
    else:
        raise NotImplementedError
    assert base_ref != base_alt
    variant["base_ref"] = base_ref
    if isinstance(base_alt, str):
        variant["base_alt"] = base_alt
    else:
        # TODO: indel
        variant["base_alt"] = ""
    return variant


def applyNovelVariant(
    backbone_seq: str,
    allele_seq: str,
    novel_variants: list[NovelVariant],
) -> str:
    """Apply the variant on allele_seq"""
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
    indent = " "
    for variant in novel_variants:
        if variant["skip"]:
            continue
        v = variant["variant"]
        pos = v.pos
        print(indent, f"Apply {v.ref}:{v.pos} {v.val} ({v.typ}) id={v.id}")
        print(indent, "  Reference", backbone_seq[pos - 5 : pos + 6])
        print(indent, "  Before:  ", allele_seq[pos - 5 : pos + 6])
        if v.typ != "single":
            print(indent, "  -> Skip (Not implement indel)")
            variant["skip"] = True
            variant["skip_reason"] = "Not implement indel"
            continue
        print(
            indent,
            "  Variant  ",
            f"     {variant['base_ref']}      -> {variant['base_alt']}",
        )

        # apply
        allele_seq = allele_seq[:pos] + variant["base_alt"] + allele_seq[pos + 1 :]
        print(indent, "  -> After:", allele_seq[pos - 5 : pos + 6])

    return allele_seq


def groupReadToBam(
    input_bam: str, output_bam: str, assign_reads: GroupPairRead
) -> None:
    """Write the assigned/grouped reads into bam file"""
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
    """Query the base on the specific position"""
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


def countFilterPileup(
    pileup_read_base: dict[str, str], reads: list[PairRead]
) -> dict[str, int]:
    """
    Filter the reads in pileup_read_base by id of given reads and
    count the occurance of each base
    """
    selected_read_id = set(read.l_sam.split("\t", 1)[0] for read in reads)
    selected_read_base = {
        id: base for id, base in pileup_read_base.items() if id in selected_read_id
    }
    return Counter(selected_read_base.values())


def splitReadsByAlleles(
    pn_typing_model: TypingWithPosNegAllele,
    predict_alleles: list[str],
) -> Iterable[tuple[str, tuple[str, ...], list[PairRead], dict[str, Variant]]]:
    """Assign the reads to alleles"""
    for gene, reads in pn_typing_model._gene_reads.items():
        typ = AlleleTyping(reads, pn_typing_model._gene_variants[gene], no_empty=False)
        assert typ.probs.shape[0] == len(reads)
        assign_reads = groupReadByAllele(typ, predict_alleles, reads)
        for alleles, reads in assign_reads.items():
            yield gene, alleles, reads, typ.variants


def discoverNovel(
    variant_name: str,
    msa_name: str,
    result_name: str,
    output_name: str,
    novel_descr: TextIO = sys.stdout,
    apply: bool = True,
) -> None:
    """
    Find novel variant (variant_name) against the allele (In result_name).
    Via pileup (bamfile derived from variant_name).

    Set apply to True to enable typing novel allele by the novel variants.
    """
    # read alleles
    result = pd.read_csv(result_name + ".tsv", sep="\t")
    predict_alleles = result["alleles"][0].split("_")
    print(predict_alleles)
    # read bam file
    bam_name = variant_name + ".no_multi"
    bamfile = AlignmentFile(bam_name + ".bam", "rb")
    data = TypingWithPosNegAllele(variant_name + ".json")
    # backbone and allele sequences
    msas: dict[str, Genemsa] = {}  # saved all genes' msa
    # init data
    allele_reads: dict[tuple[str, ...], list[PairRead]] = {}  # {allele: reads}
    allele_novel_variants: list[NovelVariant] = []  # list of novel variants
    allele_called_seqs: list[SeqRecord] = []  # new allele sequence after apply
    allele_count: dict[str, int] = defaultdict(int)
    for gene, alleles, reads, variants in splitReadsByAlleles(data, predict_alleles):
        allele_reads[alleles] = reads
        if len(alleles) > 1:
            continue
        allele = alleles[0]
        allele_count[gene] += 1
        print(f"{gene} - {allele}", file=novel_descr)

        # msa
        if gene not in msas:
            msa_gene_name = msa_name + "." + gene.split("*")[0]
            msa = Genemsa.load_msa(msa_gene_name + ".fa", msa_gene_name + ".json")
            msas[gene] = msa
        else:
            msa = msas[gene]
        allele_seq = str(msa.get(allele).replace("E", "-"))
        allele_seq_nogap = allele_seq.replace('-', '')
        backbone_seq = msa.get_reference()[1]
        print(f"  Length: {len(allele_seq_nogap)}", file=novel_descr)
        print(f"  Avg depth: {len(reads) / len(allele_seq_nogap) * 150 * 2}",
              file=novel_descr)

        # list variants
        confusion = statNovelConfusion(allele, reads, variants)
        print("  Total reads:", len(reads), file=novel_descr)
        print("  Total variants:", confusion["total"], file=novel_descr)
        for stat, c in confusion.items():
            print(f"    {stat}:", c, file=novel_descr)

        # count
        novel_variants: list[NovelVariant] = []
        stat_novel_variants = extractNovelVariant(allele, reads, variants)
        for stat, variants_count in stat_novel_variants.items():
            for variant, c in variants_count.items():
                novel_variants.append(
                    {
                        "gene": gene,
                        "allele": allele,
                        "allele_count": allele_count[gene],
                        "type": stat,
                        "variant": variant,
                        "pos": int(variant.pos),
                        "count": c,
                        "skip": False,
                        "skip_reason": "",
                        "base_ref": "",
                        "base_alt": "",
                        "pileup": {},
                    }
                )

        # filter
        threshold = 3
        for nv in novel_variants:
            if nv["count"] < threshold:
                nv["skip"] = True
                nv["skip_reason"] = "Number of variant too low"

        # get pileup from bam
        for nv in novel_variants:
            if nv["skip"]:
                continue
            pileup_read_base = queryPileup(bamfile=bamfile, ref=gene, pos=nv["pos"])
            count_pileup = countFilterPileup(pileup_read_base, reads)
            nv["pileup"] = count_pileup
            if not count_pileup:
                nv["skip"] = True
                nv["skip_reason"] = "Pileup empty"

        # filter the variant by pileup (alt depth > ref depth)
        for nv in novel_variants:
            if nv["skip"]:
                continue
            updateBaseRefAlt(nv, backbone_seq, allele_seq)
            if nv["pileup"].get(nv["base_alt"], 0) < max(nv["pileup"].values()):
                nv["skip"] = True
                nv["skip_reason"] = "ALT depths < REF depths"

        # print
        if any(map(lambda nv: not nv["skip"], novel_variants)):
            print("  List", file=novel_descr)
        for nv in novel_variants:
            if nv["skip"]:
                continue
            v = nv["variant"]
            print(
                f"    {nv['type']:5s} {v.ref}:{v.pos} {v.val} ({v.typ}) id={v.id}"
                f" num={nv['count']} bamPile={nv['pileup']}",
                file=novel_descr,
            )
        allele_novel_variants.extend(novel_variants)

        if apply:
            allele_seq = applyNovelVariant(backbone_seq, allele_seq, novel_variants)
            novel_variants = [nv for nv in novel_variants if not nv["skip"]]
            print("  Apply variant num:", len(novel_variants), file=novel_descr)

            allele_name = allele + "".join(
                f"-{nv['pos']}{nv['base_alt']}" for nv in novel_variants
            )
            allele_called_seqs.append(
                SeqRecord(
                    Seq(allele_seq),
                    id=allele_name,
                    description=",".join(
                        f"{allele}:{nv['pos']}{nv['base_ref']}>{nv['base_alt']}"
                        for nv in novel_variants
                    ),
                )
            )

    df = pd.DataFrame(allele_novel_variants)
    df["variant_type"] = [nv.typ for nv in df["variant"]]
    df["variant_id"]   = [nv.id  for nv in df["variant"]]
    df["variant_val"]  = [nv.val for nv in df["variant"]]
    df = df.drop("variant", axis=1)
    print(df, file=novel_descr)
    df.to_csv(output_name + ".variant.tsv", index=False, sep="\t")

    if apply:
        pd.DataFrame([{
            "name": output_name,
            "alleles": "_".join(i.id for i in allele_called_seqs),
        }]).to_csv(output_name + ".tsv", sep="\t", index=False)
        SeqIO.write(allele_called_seqs, output_name + ".fa", "fasta")
        groupReadToBam(bam_name + ".bam", output_name + ".bam", allele_reads)
