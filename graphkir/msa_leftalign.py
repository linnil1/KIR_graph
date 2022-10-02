"""
Left align the msa (Optional)
"""
import heapq
from typing import Iterator
from functools import reduce
from dataclasses import dataclass

from pyhlamsa import msaio, Genemsa

from .kir_msa import readFromMSAs, saveAllMsa


@dataclass(order=True)
class Segment:
    pos: int
    leng: int


def findDeletePos(seq: str) -> Iterator[Segment]:
    """
    Find the deletion position(index) and length

    Example:
      input: `AT--ATC-C`
      return: (2,2), (6,1)
    """
    pos = 0
    length = 0
    for i in range(len(seq)):
        if seq[i] == "-":
            length += 1
        else:
            if length:
                yield Segment(pos, length)
            pos = i + 1
            length = 0
    if length:
        yield Segment(pos, length)


def diff(ref_seq: str, seq: str) -> list[bool]:
    """ return two string is identical per base """
    return list(map(lambda i: i[0] == i[1], zip(ref_seq, seq)))


def findShift(ref_seq: str, seq: str, seg: Segment) -> tuple[Segment, int] | None:
    """
    Find the possible left align shift value for the deletion locus (pos + length)

    Return:
      tuple[shiftable_segment, shift]: The segment on the seg and left shift by `shift`
    """
    # slow but work
    for shift in range(seg.leng, 0, -1):
        if seg.pos - shift < 0:
            continue
        for seg_leng in range(shift, seg.leng + 1):
            # sequence at pos - shift == sequence at pos
            # or
            # diff of sequence and ref_seq at pos - shift == diff at pos
            # and
            # not overlapping
            if ("-" not in seq[seg.pos - shift : seg.pos - shift + seg_leng] and
                diff(ref_seq[seg.pos - shift : seg.pos - shift + seg_leng],
                         seq[seg.pos - shift : seg.pos - shift + seg_leng]) ==
                diff(ref_seq[seg.pos         : seg.pos         + seg_leng],
                         seq[seg.pos - shift : seg.pos - shift + seg_leng])):
                return Segment(seg.pos, seg_leng), shift

    # cannot found
    return None


def applyShift(seq: str, shift_seg: Segment, shift: int) -> str:
    """ Apply shift """
    new_seq = seq[:shift_seg.pos  - shift] \
              + seq[shift_seg.pos:         shift_seg.pos + shift_seg.leng] \
              + seq[shift_seg.pos - shift: shift_seg.pos] \
              + seq[shift_seg.pos                        + shift_seg.leng:]
    return new_seq


def leftAlign(ref_seq: str, ori_seq: str) -> str:
    """
    Left align the sequence

    Example:
      Case1:
        ```
        ACTACCACCACC
        ACTTCC---ACC
        ACT---TCCACC
        ```
      Case2:
        ```
        ACCATATATACC
        ACCATAT--ACC
        ACC--ATATACC
        ```
      Case3:
        ```
        ACCATATATACC
        ACCAT----ACC
        ACC----ATACC
        ```
      Case4:
        ```
        ACCGCCACCACC
        ACCTCC---ACC
        ACC---TCCACC
        ```
      Case5:
        ```
        ACCATATATTACC
        ACCATAT---ACC
        ACC--ATA-TACC
        ```
    """
    # deletions segments
    seg_queue = list(findDeletePos(ori_seq))
    heapq.heapify(seg_queue)

    seq = ori_seq
    while len(seg_queue):
        seg = heapq.heappop(seg_queue)
        shift_data = findShift(ref_seq, seq, seg)
        if shift_data is None:
            continue
        shift_seg, shift = shift_data
        # debug
        # print(shift_seg)
        # print(ref_seq[shift_seg.pos - shift - 5: shift_seg.pos + shift_seg.leng + 5])
        # print(    seq[shift_seg.pos - shift - 5: shift_seg.pos + shift_seg.leng + 5])
        seq = applyShift(seq, shift_seg, shift)
        # print(    seq[shift_seg.pos - shift - 5: shift_seg.pos + shift_seg.leng + 5])

        # after shift -> add into deletion queue
        heapq.heappush(seg_queue, Segment(
            shift_seg.pos - shift,
            shift_seg.leng
        ))
        # The segment is splited
        if seg.leng == shift_seg.leng:
            continue
        heapq.heappush(seg_queue, Segment(
            seg.pos + shift_seg.leng,
            seg.leng - shift_seg.leng
        ))
    # double check
    assert seq.replace('-', '') == ori_seq.replace('-', '')
    return seq


def msaLeftAlign(msa_ori: Genemsa) -> Genemsa:
    """ Left align the MSA (block-wise) """
    msa_aligned = []
    for msa in msa_ori.split():
        for name, seq in list(msa.alleles.items()):
            ref_seq = msa.get_reference()[1]
            msa.alleles[name] = leftAlign(ref_seq, seq)
            # debug
            # print(name)
            # msa.alleles[name + "_a"] = leftAlign(ref_seq, seq)
        # msa = msa.sort_name()
        # print(msa.format_alignment_diff())
        msa_aligned.append(msa)
    return reduce(lambda i, j: i + j, msa_aligned)


def genemsaLeftAlign(input_prefix: str, output_prefix: str):
    """
    Left align our MSA in `input_prefix.{gene}` and save in `output_prefix.{gene}`
    Also reconstruct our backbone sequence
    """
    msas = readFromMSAs(input_prefix)
    new_msas = {}
    for gene, msa in msas.items():
        print("left_align", gene)
        refname = msa.get_reference()[0]
        assert refname == f"{gene}*BACKBONE"
        msa = msaLeftAlign(msa)
        msa = msa.remove(refname)
        new_msas[gene] = msa
    saveAllMsa(new_msas, output_prefix)


if __name__ == "__main__":
    # test
    cases = [
        ("ACTACCACCACC",
         "ACTTCC---ACC",
         "ACT---TCCACC"),
        ("ACCATATATACC",
         "ACCATAT--ACC",
         "ACC--ATATACC"),
        ("ACCATATATACC",
         "ACCAT----ACC",
         "ACC----ATACC"),
        ("ACCGCCACCACC",
         "ACCTCC---ACC",
         "ACC---TCCACC"),
        ("ACCATATATTACC",
         "ACCATAT---ACC",
         "ACC--ATA-TACC"),
    ]

    for ref, alt, expect in cases:
        print(ref)
        print(alt)
        print(leftAlign(ref, alt))
        assert leftAlign(ref, alt) == expect
        print()
    """
    gene = "KIR2DL1S1"
    genemsaLeftAlign(f"index5/kir_2100_ab_2dl1s1.{gene}",
                     f"index5/kir_2100_ab_2dl1s1.leftalign.{gene}")
    """
