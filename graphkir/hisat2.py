"""
index + bam -> JSON

* backbone's read
* read's variant (postive, negative)
* variant's allele
"""
import re
import copy
import json
import bisect
from typing import TypedDict, Iterable
from dataclasses import dataclass, field, asdict

from .utils import runDocker, samtobam
from .msa2hisat import Variant
from .pileup import PileupCount, getPileupBaseRatio

# TODO
# = typing_common.identify_ambigious_diffs(ref_seq,


@dataclass
class PairRead:
    """
    Save the record and its variants of each variant in pair read

    Attributes:
      l_sam: The read records
      r_sam: The read records which is the other read in the pair
      multiple: Number of mapping can this pair mapped on references
      backbone: reference name
      lpv: positive variants contains in the l_sam
      rpv: positive variants contains in the r_sam
      lnv: negative variants contains in the l_sam
      rnv: negative variants contains in the r_sam
    """

    # records strings
    l_sam: str = ""
    r_sam: str = ""

    # some usful info
    multiple: int = 1
    backbone: str = ""

    # left/right's positive/negative variants
    lpv: list[str] = field(default_factory=list)
    lnv: list[str] = field(default_factory=list)
    rpv: list[str] = field(default_factory=list)
    rnv: list[str] = field(default_factory=list)


class ReadsAndVariantsData(TypedDict):
    """
    Result data after parsing bam_file

    Attributes:
      variants: variants data in the graph (including novel variants of the sample)
      reads: The positive/negative variants of each pair read
    """

    variants: list[Variant]
    reads: list[PairRead]


def hisatMap(index: str, f1: str, f2: str, output_file: str, threads: int = 1) -> None:
    """run hisat2"""
    assert output_file.endswith(".bam")
    output_name = output_file.rsplit(".", 1)[0]
    runDocker(
        "hisat",
        f"""\
        hisat2 --threads {threads} -x {index} -1 {f1} -2 {f2} \
               --no-spliced-alignment --max-altstried 64 --haplotype \
               -S {output_name}.sam
        """,
    )
    samtobam(output_name)


def getNH(sam_info: str) -> int:
    """Extract NH from record"""
    NH_re = re.findall(r"NH:i:(\d+)", sam_info)
    if NH_re:
        return int(NH_re[0])
    return 1


def readBam(bam_file: str) -> list[str]:
    """Read bam file via samtools"""
    proc = runDocker(
        "samtools", f"samtools sort -n {bam_file} -O SAM", capture_output=True
    )
    return str(proc.stdout).split("\n")


def readBamHeader(bam_file: str) -> str:
    """Read header of bam file via samtools"""
    proc = runDocker("samtools", f"samtools view -H {bam_file}", capture_output=True)
    return str(proc.stdout)


def readLink(index: str) -> dict[str, list[str]]:
    """
    Read hisat2 .link file

    Returns:
        Dictionary[variant_id, list of allele related to the variant]
    """
    id_alleles = {}
    with open(index + ".link") as f:
        for line in f:
            # hv25	KIR2DL2*0010101 KIR2DL2*0010103
            allele_id, alleles = line.strip().split("\t")
            id_alleles[allele_id] = alleles.split()
    return id_alleles


def readExons(index: str) -> dict[str, list[tuple[int, int]]]:
    """
    Read exon's information via hisat2 .locus file

    Returns:
        Dictionary[reference, list[exon_start_position, exon_end_position]]
    """
    # locus = 1-base position
    # variant = 0-base position
    # KIR2DL1*BACKBONE        KIR2DL1*BACKBONE        0       14761   14761
    # 269-303 1267-1303 2031-2313 3754-4054 5588-5882 +
    gene_exons = {}
    with open(index + ".locus") as f:
        for line in f:
            gene, _, _, _, _, exons_str, _ = line.split("\t")
            gene_exons[gene] = [
                (int(i.split("-")[0]) - 1, int(i.split("-")[1]) - 1)
                for i in exons_str.split(" ")
            ]
    return gene_exons


def readVariants(index: str) -> list[Variant]:
    """
    Read hisat2 .snp file

    Returns:
        A list of variants (only basic information)
    """
    variants = []
    with open(index + ".snp") as f:
        for line in f:
            # hv7	single	KIR2DL2*BACKBONE	5819	T
            # hv8	deletion	KIR2DL2*BACKBONE	6626	2
            variant_id, var_type, ref, pos, val = line.strip().split("\t")
            v = Variant(
                typ=var_type,
                ref=ref,
                pos=int(pos),
                id=variant_id,
                val=val if var_type != "deletion" else int(val),
            )
            variants.append(v)
    return variants


def getVariants(index: str) -> list[Variant]:
    """
    Read hisat2 files to reconstruct variants information

    Returns:
        A list of variants (all information)
    """
    # basic info
    variants = readVariants(index)

    # list of alleles related to variant
    id_alleles = readLink(index)
    for v in variants:
        assert v.id
        v.allele = id_alleles.get(v.id, [])

    # if exon in pos
    exon_region = readExons(index)
    for v in variants:
        v.in_exon = isInExon(exon_region[v.ref], v)
    return sorted(variants)


def isInExon(exons: list[tuple[int, int]], variant: Variant) -> bool:
    """
    Is the variant in exon region

    Args:
        exons: The list of exons start/end position from `readExons`
    Returns:
        retrun True if the variant in exon
    """
    # TODO: binary search, but naive method is not slow
    for exon in exons:
        if exon[0] <= variant.pos < exon[1]:
            return True
        if (
            variant.typ == "deletion"
            and variant.pos < exon[0]
            and variant.pos + variant.val >= exon[0]  # type: ignore
        ):
            return True
    return False


def readPair(bam_file: str) -> Iterable[tuple[str, str]]:
    """
    Extract paired read from bam

    Args:
      bam_file: The bam file path

    Yields:
      left_record(str), right_record(str)
    """
    num_reads = 0
    num_pairs = 0
    reads = {}  # type: dict[tuple[str, str, str], str]  # A temporary dict

    for line in readBam(bam_file):
        # skip header
        if not line or line.startswith("@"):
            continue
        # skip accidentally pipe stderr to stdout
        if line.startswith("[bam_sort_core]"):
            continue

        # preprocess
        read_id, flag_str, ref, pos, _, _, next_ref, next_pos = line.split("\t")[:8]
        if next_ref != "=":
            continue
        num_reads += 1
        flag = int(flag_str)

        # If pair read in temp dict -> remote from temp, check and yield
        next_ref = ref
        this_id = (read_id, ref,      pos,      flag & 256)
        next_id = (read_id, next_ref, next_pos, flag & 256)
        if next_id in reads:
            next_line = reads[next_id]
            # is left and right
            if ((int(next_line.split("\t")[1]) | flag) & (64 + 128)) != 64 + 128:
                print("Strange case", line, next_line)
                continue
            # success
            del reads[next_id]
            num_pairs += 1
            yield line, next_line
        # if not find -> save in temp
        else:
            reads[this_id] = line

    # summary
    print("Reads:", num_reads, "Pairs:", num_pairs)


def recordToRawVariant(line: str) -> tuple[list[Variant], list[int]]:
    """
    Turn bam record into a list of variant (ref + pos + alt)

    We use these bam information

    * pos (Start position)
    * cigar (diff from reference)

        Example: `98M1D52M`

        Format: [(cigar_length, cigar_type)]

    * MD (info when walking on graph)

        Example: `MD:Z:43^G19^TGGAGATATGGGCCTGGGTG88`

        Format:

        * if number: match number
        * elif starwith '^":  deletion (and it's deleted bases)
        * else mismatch
        * Insertion is not included in this format

    * Zs (Annotation of the variant)

        Example: `Zs:Z:43|D|hv862,19|D|hv868`

        Format:

        * (gap|type|var_id)
        * gap shows the gap between each var_id
        * if type is S (Single), the var_id has width = 1, otherwise = 0

    When implement, we need to maintain these five variable in the same time

    * start, pos:                Absolute Position to reference
    * read_i:                    Position of sequence of read
    * zs_i, zs_pos:              The Zs index and position
    * md_i, md_len, md_len_last: The md index and position
    * cigar:                     Cigar is mainly couple with md

    And here is their relationship

    ```
    pos = start + cigar_length * n
        = start + read_i
        = start + zs_gap * n + 1 * (if single)

    # When iterate the md
    md_pos      = pos + md_len_last
    current_pos = pos + md_len
    # (md_len will always < cigar_length, otherwise pos will increase)
    # and md_len > md_len_last i.e. current_pos > md_pos
    ```

    Args:
        line: bam record
    Returns:
        * The list of variant
        * softclip tuple[head clipping length, tail clipping length]
    """
    # read info
    cols       = line.strip().split('\t')
    backbone   = cols[2]
    start      = int(cols[3]) - 1  # no shift from locus
    cigar_strs = re.findall(r'\d+\w', cols[5])
    cigars     = map(lambda i: (i[-1], int(i[:-1])), cigar_strs)
    read_seq   = cols[9]
    # read_qual  = cols[10]
    zs         = readZs(cols[11:])
    md         = readMd(cols[11:])
    # return print(list(cigars), list(zs), list(md))

    # for inner function
    pos        = start  # The absolute position (default=start)
    read_i     = 0      # Read position (relative to start)
    zs_i, md_i = 0, 0   # Array index for zs and md
    zs_pos     = 0      # Accumulated zs (relative to start)
    md_len     = 0      # The len of this md to pos

    def findZs(var_type: str) -> str:
        nonlocal zs_i, zs_pos
        var_id = "unknown"
        if (
            zs_i < len(zs)
            and read_i + md_len == zs_pos + zs[zs_i][0]
            and zs[zs_i][1] == var_type
        ):
            zs_pos += zs[zs_i][0]
            if var_type == "S":
                zs_pos += 1
            var_id = zs[zs_i][2]
            zs_i += 1
        return var_id

    def mdMatch() -> Iterable[Variant]:
        """
        Read md(s) for the cigar

        * Case1: md_len < cigar_length
                 cigar: 150M
                 md: 77 A 1 C 70
        * Case2: md_len > cigar_length
                 cigar: 1I149M
                 md: 19 A 130
        """
        nonlocal md_len, md_i
        md_len_last = 0  # The length from last md to pos
        while True:
            # print("Match", md_i, md_len, md_len_last)
            # read md
            if md_len <= md_len_last and md_i < len(md) and type(md[md_i]) is int:
                md_len += int(md[md_i])
                md_i += 1

            # md_len > cigar
            # leave md_len !=0 for latter used
            # the sequence before current_pos is the same (after md_len_last)
            if md_len >= cigar_length:
                md_len -= cigar_length
                yield Variant(
                    typ="match",
                    ref=backbone,
                    pos=pos + md_len_last,
                    length=cigar_length - md_len_last,
                )
                break

            # read the base
            # * md[md_i] is reference
            # * read_seq[pos] is alt

            read_base = read_seq[read_i + md_len]
            # why cannot remove 0 before: Think about this case MD:Z:7A37A0G55T4^CA0T42
            if md[md_i] == 0:  # remove 0
                md_i += 1
            assert str(md[md_i]) in "ACGT"
            assert str(md[md_i]) != read_base
            md_i += 1

            # the sequence before current_pos is the same
            if md_len > md_len_last:
                yield Variant(
                    typ="match",
                    ref=backbone,
                    pos=pos + md_len_last,
                    length=md_len - md_len_last,
                )

            # the current position is mismatch
            yield Variant(
                typ="single",
                ref=backbone,
                pos=pos + md_len,
                length=1,
                val=read_base,
                id=findZs("S"),
            )
            md_len += 1
            md_len_last = md_len

            # cigar OK, increase the pos and set md_len = md_len_last = 0
            if md_len == cigar_length:
                md_len = 0
                break

    def mdInsertion() -> list[Variant]:
        # md does not contains insertion info
        return [
            Variant(
                typ="insertion",
                ref=backbone,
                pos=pos,
                val=read_seq[read_i : read_i + cigar_length],
                length=cigar_length,
                id=findZs("I"),
            )
        ]

    def mdDeletion() -> list[Variant]:
        nonlocal md_i
        assert md[md_i] == "^"
        md_i += 1
        while md_i < len(md) and type(md[md_i]) is not int and str(md[md_i]) in "ACGT":
            md_i += 1
        return [
            Variant(
                typ="deletion",
                ref=backbone,
                pos=pos,
                val=cigar_length,
                length=cigar_length,
                id=findZs("D"),
            )
        ]

    cmp_list: list[Variant] = []
    soft_clip = [0, 0]
    # print(line)
    for cigar_i, (cigar_op, cigar_length) in enumerate(cigars):
        # print("cigar", cigar_op, cigar_length)
        # print("pos read_i", pos, read_i)
        # print("ZS", zs, zs_i, zs_pos)
        # print('MD', md, md_i, md_len)
        # print('cmp_list', cmp_list)
        if md_i < len(md) and md[md_i] == 0:  # remove 0
            md_i += 1
        if cigar_op == "M":
            cmp_list.extend(mdMatch())
        elif cigar_op == "I":
            cmp_list.extend(mdInsertion())
        elif cigar_op == "D":
            cmp_list.extend(mdDeletion())
        elif cigar_op == "S":
            if cigar_i == 0:
                soft_clip[0] = cigar_length
            else:
                soft_clip[1] = cigar_length
            zs_pos += cigar_length
        elif cigar_op == "N":
            raise NotImplementedError("Cannot typing with splicing")
        else:
            raise NotImplementedError

        if cigar_op in "MND":
            pos += cigar_length
        if cigar_op in "MIS":
            read_i += cigar_length

    if md_i < len(md) and md[md_i] == 0:  # remove 0
        md_i += 1
    # double check
    assert zs_i == len(zs)
    assert md_i == len(md)
    assert read_i == len(read_seq)
    return cmp_list, soft_clip


def readZs(cols: list[str]) -> list[tuple[int, str, str]]:
    """Read hisat2 Zs infomation"""
    for col in cols:
        if col.startswith("Zs"):
            # Zs:Z:43|D|hv862,19|D|hv868
            zs = col[5:].split(",")
            zs1 = map(lambda i: i.split("|"), zs)
            zs2 = map(lambda i: (int(i[0]), i[1], i[2]), zs1)
            return list(zs2)
    return []


def readMd(cols: list[str]) -> list[int | str]:
    """Read hisat2 MD infomation"""
    for col in cols:
        if col.startswith("MD"):
            # MD:Z:43^G19^TGGAGATATGGGCCTGGGTG88
            md = re.findall(r"(\d+|.)", col[5:])
            md1 = map(lambda i: int(i) if i.isdigit() else i, md)  # type: ignore
            return list(md1)
    return []


def filterRead(line: str, num_editdist: int = 4) -> bool:
    """
    Filter reads (same in hisat2)

    * Concordently mapped
    * edit distance < num_editdist

    Returns:
        If return True, the read can pass the criteria
    """
    flag = int(line.split("\t")[1])

    # Concordantly mapped?
    if flag & 2 == 0:
        # print("Not concordantly mapped read")
        # print(line)
        return False

    # unmapped(not needed)
    # if flag & 4 or flag & 8:
    #     return False

    # read additional flag
    NM = None
    for col in line.strip().split("\t")[11:]:
        if col.startswith("NM"):
            NM = int(col[5:])

    if NM is None or NM > num_editdist:
        # print(f"Read > editdist (NM = {NM})")
        # print(line)
        return False

    # secondary
    # if flag & 256:
    #     print("Secondary Read")
    #     return False
    return True


def findVariantId(variant: Variant, variants_map: dict[Variant, Variant]) -> Variant:
    """
    Try to match the known variant in hisat2

    Maybe the variant is not located in .snp.index but in .snp

    If this is a novel variant, the variant will add into variants_map

    Args:
      variant: The variant for searching
      variants_map: All variants in HISAT2
    """
    if variant in variants_map:
        return variants_map[variant]

    # cannot find: set it novel
    if variant.typ in ["single", "insertion", "deletion"]:
        # print("Cannot find", variant)
        variant.id = f"nv{Variant.novel_id}"
        Variant.novel_id += 1
        variants_map[variant] = variant
        return variant

    # print("Other", variant)
    assert variant.typ == "match"
    return variant


def errorCorrection(variant: Variant, pileup: PileupCount) -> Variant:
    """
    Error correction of the base based on pileup

    Args:
      variant: The target variant
      pileup: Pileup information (created by `getPileupBaseRatio`)
    """
    # Note that hisat2 error correction has criteria
    # `if num_nt >= 20 and (count >= num_nt * 0.2 or count >= 7)`
    # TODO: split the match when some base in match is sequencing error
    if variant.typ != "single":
        return variant

    # SNP
    p = pileup.get((variant.ref, variant.pos))
    if not p:
        # print("error pileup", variant)
        return variant

    # insufficient depth
    if p["all"] < 20:
        return variant

    # not minority
    if p.get(str(variant.val), 0) > 0.2:
        # print("reserve", variant, p)
        return variant

    # minority -> change to majority
    # case: AAAAAAAATAAAA where REF = G
    if any(i[1] >= 0.8 for i in p.items() if i[0] != "all"):
        variant.val = max([i for i in p.items() if i[0] != "all"], key=lambda i: i[1])[0]
        return variant

    # if variant.id.startswith("hv"):
    #     print("to_match", variant, p)
    # case: AAAAAACCCCCCCCCCC REF = C
    # neither set to A or C is wrong
    # if I set it to match or novel -> it will become negative
    # variant.typ = "match"

    # minority but cannot assign to any variant
    # AAAAAAATTTTTTTC
    variant.val = "N"
    return variant


def recordToVariants(
    record: str,
    variants_map: dict[Variant, Variant],
    pileup: None | PileupCount = None,
    ignore_softclip: bool = False,
) -> list[Variant]:
    """
    Extract variants (with its id) in the read

    Args:
        record: bam record
        variants_map:
            All the variants in the graph.

            Note that if there is novel variant, the variant will add into variants_map
        pileup:
            created by `getPileupBaseRatio`,
            if not provided, it'll not execute error correction
        ignore_softclip:
            Set it False, it'll treat soft clip as fail reads (return empty list)

    Returns:
        The sorted list of (positive) Variant in the read
    """
    read_variants, soft_clip = recordToRawVariant(record)
    if not ignore_softclip and sum(soft_clip) > 0:
        return []

    if pileup:
        read_variants = [errorCorrection(v, pileup) for v in read_variants]

    read_variants = [findVariantId(v, variants_map) for v in read_variants]
    return sorted(read_variants)


def getVariantsBoundary(
    read_variants: list[Variant], variants: list[Variant]
) -> tuple[int, int]:
    """
    Find the lower and upper bound of read

    i.e. find all the variants inside the read region (include negative)

    Args:
        read_variants: The variants inside the reads
        variants: All the variants in the graph

    Returns:
        [left, right) (the index of .snp variants)
    """
    left = read_variants[0].pos
    rigt = read_variants[-1].pos + read_variants[-1].length
    ref = read_variants[0].ref
    return (
        bisect.bisect_left(variants, Variant(ref=ref, pos=left, typ="single", val="A")),
        bisect.bisect_left(variants, Variant(ref=ref, pos=rigt, typ="single", val="T")),
    )


def getPNFromVariantList(
    read_variants: list[Variant],
    variants: list[Variant],
    exon_only: bool = False,
    discard_novel_index: bool = True,
) -> tuple[list[Variant], list[Variant]]:
    """
    Extract positive and negative variants of read

    Note:
      * Remove novel insertion and deletion
      * Find negative variant (exclude deletion at the front and end)

    Args:
        read_variants: The **sorted** variants in the read
        variants: The **sorted** variants list of the hisat2
        exon_only: Extract variants only on exon
        discard_novel_index: set it True if assume novel insert/deletion = mapping error

    Returns:
        * list of positive allele
        * list of negative allele
    """
    # Find all variant in the region
    if not read_variants:
        return [], []
    left, right = getVariantsBoundary(read_variants, variants)
    pos_right = read_variants[-1].pos + read_variants[-1].length
    assert left <= right

    if discard_novel_index:
        if any(
            v.typ == "insertion" and v.id.startswith("nv")  # type: ignore
            for v in read_variants
        ):
            return [], []
        if any(
            v.typ == "deletion" and v.id.startswith("nv")  # type: ignore
            for v in read_variants
        ):
            return [], []

    # exclude for negative
    exclude_variants = set()

    # exclude cannot error correction variant
    for v in read_variants:
        if v.val == "N":
            for i in "ATCG":
                va = copy.deepcopy(v)
                va.val = i
                exclude_variants.add(va)

    # remove match
    read_variants = [v for v in read_variants if v.typ != "match"]
    # remove novel
    # read_variants = [v for v in read_variants if v.id and not v.id.startswith("nv")]

    # positive variation
    if exon_only:
        positive_variants = [v for v in read_variants if v.in_exon]
    else:
        positive_variants = read_variants
    exclude_variants.update(positive_variants)

    # print('positive', read_variants)

    # negative
    negative_variants = []
    for v in variants[left:right]:
        if v in exclude_variants:
            continue
        if exon_only and not v.in_exon:
            continue
        # Deletion should not over the right end of read
        # Because some deletion is ambiguous until last one
        # TODO: change 10 -> to repeat length
        if v.typ == "deletion" and v.pos + v.val + 10 >= pos_right:  # type: ignore
            continue
        # print('negative', v)
        negative_variants.append(v)
        # TODO: maybe some deletion is before left
        # they use gene_var_maxrights to deal with

    return positive_variants, negative_variants


def extractVariant(
    pair_reads: Iterable[tuple[str, str]],
    variants: list[Variant],
    pileup: None | PileupCount = None,
) -> ReadsAndVariantsData:
    """bam records (pair) -> variants information per read"""
    # this dictionary may changes
    variants_map = {v: v for v in variants}

    # init
    reads = []
    tot_pair = 0
    for left_record, right_record in pair_reads:
        tot_pair += 1
        # main: extract bam -> variant
        # main: + annotate ID
        lv = recordToVariants(left_record, variants_map, pileup)
        rv = recordToVariants(right_record, variants_map, pileup)

        # (left, right) x (positive, negative)
        lp, ln = getPNFromVariantList(lv, variants)
        rp, rn = getPNFromVariantList(rv, variants)

        # save the allele information for the read
        reads.append(
            PairRead(
                lpv=[v.id for v in lp if v.id is not None],
                lnv=[v.id for v in ln if v.id is not None],
                rpv=[v.id for v in rp if v.id is not None],
                rnv=[v.id for v in rn if v.id is not None],
                l_sam=left_record,
                r_sam=right_record,
                multiple=getNH(left_record),
                backbone=left_record.split("\t")[2],
            )
        )
    print("Filterd pairs", tot_pair)
    reads_data: ReadsAndVariantsData = {
        "variants": list(variants_map.values()),
        "reads": reads,
    }
    return reads_data


def writeReadsAndVariantsData(reads_data: ReadsAndVariantsData, filename: str) -> None:
    """Write result to json"""
    with open(filename, "w") as f:
        json.dump(
            {
                "variants": [asdict(i) for i in reads_data["variants"]],
                "reads": [asdict(i) for i in reads_data["reads"]],
            },
            f,
        )


def loadReadsAndVariantsData(filename: str) -> ReadsAndVariantsData:
    """Load result from json"""
    with open(filename) as f:
        reads_data = json.load(f)
    return {
        "variants": [Variant(**i) for i in reads_data["variants"]],
        "reads": [PairRead(**i) for i in reads_data["reads"]],
    }


def saveSam(filename: str, header: str, reads: Iterable[PairRead]) -> None:
    """
    save records into sam file

    Args:
      filename: samfile filename
      header: samfile header
      reads: samfile records
    """
    with open(filename, "w") as f:
        f.writelines(header)
        f.writelines(map(lambda i: i.l_sam + "\n" + i.r_sam + "\n", reads))


def saveReadsToBam(
    reads_data: ReadsAndVariantsData,
    filename_prefix: str,
    bam_file: str,
    filter_multi_mapped: bool = False,
) -> None:
    """
    Save the reads into sam/bamfile (`{filename_prefix}.bam`)

    Args:
        bam_file: The header of samfile is from the old bamfile `bam_file`
        filter_multi_mapped: Remove multiple mapped reads
    """
    print(f"Save bam in {filename_prefix}.bam")
    sam_header = readBamHeader(bam_file)
    reads = reads_data["reads"]  # type: Iterable[PairRead]
    if filter_multi_mapped:
        reads = filter(lambda r: r.multiple == 1, reads_data["reads"])
    saveSam(filename_prefix + ".sam", sam_header, reads)
    samtobam(filename_prefix)


def extractVariantFromBam(
    index: str, bam_file: str, output_prefix: str, error_correction: bool = True
) -> None:
    """
    Extract reads and variants from bamfile

    1. Filter bad reads
    2. Call the variant from the read
    3. Annotate the variant by known variants (in index)
    4. Save variants and reads information into json
    5. Save reads into bam again
    6. Save multiple-mapped-reads-removed bamfile

    Args:
        bam_file: The bamfile for extraction
        index: The hisat2 format (no `.graph` if suffix)
        output_prefix: The prefix of output filename
    """
    # read file
    variants = getVariants(index)
    pair_reads = readPair(bam_file)
    pair_reads = filter(lambda lr: filterRead(lr[0]) and filterRead(lr[1]), pair_reads)
    if error_correction:
        pileup = getPileupBaseRatio(bam_file)
    else:
        pileup = None

    # main
    reads_data = extractVariant(pair_reads, variants, pileup=pileup)
    print(f"Save allele per reads in {output_prefix}.json")
    writeReadsAndVariantsData(reads_data, f"{output_prefix}.json")

    # save to another format
    saveReadsToBam(reads_data, output_prefix              , bam_file)
    saveReadsToBam(reads_data, output_prefix + ".no_multi", bam_file, filter_multi_mapped=True)


def removeMultipleMapped(reads_data: ReadsAndVariantsData) -> ReadsAndVariantsData:
    """actually remove NH != 1 reads"""
    return {
        "variants": reads_data["variants"],
        "reads": list(filter(lambda i: i.multiple == 1, reads_data["reads"])),
    }
