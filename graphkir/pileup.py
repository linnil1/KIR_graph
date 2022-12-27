"""
Pileup utilies
"""
import re
from collections import Counter
from typing import Iterator

from .utils import runDocker

PileupCount = dict[tuple[str, int], dict[str, float]]


def parsePileupBase(bases: str) -> Iterator[str]:
    """Extract base from mpileup fields[4]"""
    # chrM    1       N       10      ^]G^]G^KG^OG^]G^]G^TG^(G^KG^VG  DDDShDDDDD
    # chr1    17344300        N       2       A$^]A   <=
    # strange case: ??
    # chr1    17351906        N       1       *       E
    i = 0
    while i < len(bases):
        if bases[i] == "$":  # end of seq
            i += 1
            continue
        if bases[i] == "*":  # deletion
            yield bases[i]
            i += 1
            continue
        if bases[i] in "+-":  # insertion
            a = re.findall(r"(\d+)", bases[i + 1 :])
            assert a
            i += 1 + len(a[0]) + int(a[0])
            continue
        if bases[i] == "^":  # mapping quality
            i += 2
            continue
        yield bases[i]
        i += 1


def readPileup(bam_file: str) -> Iterator[tuple[str, int, int, str]]:
    """
    Read mpileup of the bam file

    Returns:
      (reference, position, depth, bases)
    """
    proc = runDocker("samtools", f"samtools mpileup -a {bam_file}", capture_output=True)
    for line in proc.stdout.split("\n"):
        if not line or "[mpileup]" in line:
            continue
        fields = line.split("\t")
        yield fields[0], int(fields[1]) - 1, int(fields[3]), fields[4]


def getPileupBaseRatio(bam_file: str) -> PileupCount:
    """
    Calculate the base ratio in pileup on the ref + pos

    The base include {'A', 'C', 'T', 'all', '*', 'G'}
    where `all` is the depth (total count of all bases)
    and `*` is for deletion

    Returns:
        `{(ref, pos): {"A": 0.2, "C": 0.8, "all": 30}}`
    """

    stat = {}
    for ref, pos, depth, pileup_bases in readPileup(bam_file):
        if depth == 0:
            continue

        bases = "".join(parsePileupBase(pileup_bases))
        assert depth == len(bases)

        count = Counter(bases.upper())
        s = sum(count.values())
        stat[(ref, pos)] = {k: v / s for k, v in count.items()}
        stat[(ref, pos)]["all"] = s
    return stat
