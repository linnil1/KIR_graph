# This script is copy from hisat-genotype/hisatgenotype_typing_core
import re
import os
import sys
import math
import json
import bisect
import subprocess
from typing import Union, ClassVar
from collections import defaultdict
from dataclasses import dataclass, field
import numpy as np
from Bio import SeqIO

# sam
num_editdist = 5
# sam
NH_filtering = False
# em
diff_threshold = 0.0001
# em
iter_max = 100
# em
norm_by_length = False


@dataclass
class Variant:
    pos: int
    typ: str
    ref: str = None
    val: Union[int, str] = None
    allele: list = field(default_factory=list)
    length: int = 1
    in_exon: bool = False
    id: str = None

    # allele: list = field(default_factory=list)
    # freq: float = None
    order_type: ClassVar[dict] = {"insertion": 0, "single": 1, "deletion": 2}
    order_nuc: ClassVar[dict] = {"A": 0, "C": 1, "G": 2, "T": 3}
    novel_count: ClassVar[int] = 0

    def __lt__(self, othr):
        # The sorted sta
        return (self.ref, self.pos) < (othr.ref, othr.pos)

    def __eq__(self, othr):
        return (self.pos, self.ref, self.typ, self.val) \
               == (othr.pos, othr.ref, othr.typ, othr.val)

    def __hash__(self):
        return hash((self.pos, self.ref, self.typ, self.val))


def readBam(alignment_fname):
    """ get bam file via samtools in command line """
    alignview_cmd = ["samtools", "view", "-h", alignment_fname,
                     "|", "samtools", "sort", "-n", "-", "-O", "SAM"]
    proc = subprocess.Popen(" ".join(alignview_cmd),
                            shell=True,
                            universal_newlines=True,
                            stdout=subprocess.PIPE,
                            stderr=open("/dev/null", 'w'))
    return proc.stdout


def readPair(alignment_fname):
    """ Extract paired read from bam """
    num_reads = 0
    num_pairs = 0
    reads = {}  # A temporary dict

    for line in readBam(alignment_fname):
        # skip header
        if line.startswith("@"):
            continue

        # preocess
        read_id, flag, ref, pos, _, _, next_ref, next_pos = line.split('\t')[:8]
        if next_ref != "=":
            print("Skip")
            # print(line)
            continue

        # If pair read in temp dict -> remote from temp, check and yield
        num_reads += 1
        next_ref = ref
        next_id = (read_id, next_ref, next_pos)
        if next_id in reads:
            next_line = reads[next_id]
            # is left and right
            if ((int(next_line.split('\t')[1]) | int(flag)) & (64 + 128)) != 64 + 128:
                print("Strange case", line, next_line)
                continue
            # success
            del reads[next_id]
            num_pairs += 1
            yield line, next_line
        # if not find -> save in temp
        else:
            reads[(read_id, ref, pos)] = line

    # summary
    print("Reads:", num_reads, "Pairs:", num_pairs)


def filterRead(line):
    """ Filter some reads (same in hisat2) """
    flag = int(line.split('\t')[1])

    # Concordantly mapped?
    if flag & 2 == 0:
        # print("Not concordantly mapped read")
        # print(line)
        return False

    # unmapped(not needed)
    # if flag & 4 or flag & 8:
    #     return False

    # read additional flag
    NH, NM = None, None
    for col in line.strip().split('\t')[11:]:
        if col.startswith("NM"):
            NM = int(col[5:])
        elif col.startswith("NH"):
            NH = int(col[5:])

    # filter
    if NH_filtering and NH > 1:
        # print("Read not unique aligned")
        return False
    if NM > num_editdist:
        # print(f"Read > editdist (NM = {NM})")
        # print(line)
        return False

    # secondary
    # if flag & 256:
    #     print("Secondary Read")
    #     return False
    return True


def filterPair(i):
    line_l, line_r = i
    return filterRead(line_l) and filterRead(line_r)


def readZsMd(cols):
    """ Read hisat2 Zs MD infomation """
    zs, md = [], []
    for col in cols:
        if col.startswith("Zs"):
            # Zs:Z:43|D|hv862,19|D|hv868
            zs = col[5:].split(',')
            zs = map(lambda i: i.split('|'), zs)
            zs = map(lambda i: (int(i[0]), *i[1:]), zs)
        elif col.startswith("MD"):
            # MD:Z:43^G19^TGGAGATATGGGCCTGGGTG88
            md = re.findall(r'(\d+|.)', col[5:])
            md = map(lambda i: int(i) if i.isdigit() else i, md)

    return list(zs), list(md)


def record2Variant(line):
    """
    Turn bam record into a list of variant format

    Variant format = (type, pos, length, varid, seq_in_that_position)
    * type: match, single(mismatch), insertion, deletion
    * pos: int
    * length: int
    * varid: hvxxx if exist else unknown

    We use these bam information to transfer the type
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

    When implement, we need to maintain these five category in the same time
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

    # when iter the md
    # (md_len will always < cigar_length, otherwise pos will increase)
    # and md_len > md_len_last i.e. current_pos > md_pos
    md_pos      = pos + md_len_last
    current_pos = pos + md_len

    ```
    """
    # read info
    cols = line.strip().split('\t')
    start      = int(cols[3]) - 1  # no shift from locus
    backbone   = cols[2]
    read_seq   = cols[9]
    read_qual  = cols[10]
    cigars     = re.findall(r'\d+\w', cols[5])
    cigars     = map(lambda i: (i[-1], int(i[:-1])), cigars)
    zs, md     = readZsMd(cols[11:])
    # return print(list(cigars), list(zs), list(md))

    # for inner function
    pos        = start  # The absolute position (default=start)
    read_i     = 0      # Read position (relative to start)
    zs_i, md_i = 0, 0   # Array index for zs and md
    zs_pos     = 0      # Accumulated zs (relative to start)
    md_len     = 0      # The len of this md to pos

    def findZs(var_type):
        nonlocal zs_i, zs_pos
        var_id = "unknown"
        if zs_i < len(zs) and read_i + md_len == zs_pos + zs[zs_i][0] and zs[zs_i][1] == var_type:
            zs_pos += zs[zs_i][0]
            if var_type == "S":
                zs_pos += 1
            var_id = zs[zs_i][2]
            zs_i += 1
        return var_id

    def mdMatch():
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
                md_len += md[md_i]
                md_i += 1

            # md_len > cigar
            # leave md_len !=0 for latter used
            # the sequence before current_pos is the same (after md_len_last)
            if md_len >= cigar_length:
                md_len -= cigar_length
                yield Variant(typ="match",
                              ref=backbone,
                              pos=pos + md_len_last,
                              length=cigar_length - md_len_last)
                break

            # read the base
            # * md[md_i] is reference
            # * read_seq[pos] is alt

            read_base = read_seq[read_i + md_len]
            # why cannot remove 0 before: Think about this case MD:Z:7A37A0G55T4^CA0T42
            if md[md_i] == 0:  # remove 0
                md_i += 1
            assert md[md_i] in "ACGT"
            assert md[md_i] != read_base
            md_i += 1

            # the sequence before current_pos is the same
            if md_len > md_len_last:
                yield Variant(typ="match",
                              ref=backbone,
                              pos=pos + md_len_last,
                              length=md_len - md_len_last)

            # the current position is mismatch
            yield Variant(typ="single",
                          ref=backbone,
                          pos=pos + md_len,
                          length=1,
                          val=read_base,
                          id=findZs("S"))
            md_len += 1
            md_len_last = md_len

            # cigar OK, increase the pos and set md_len = md_len_last = 0
            if md_len == cigar_length:
                md_len = 0
                break

    def mdInsertion():
        # md does not contains insertion info
        return [Variant(typ="insertion",
                        ref=backbone,
                        pos=pos,
                        val=read_seq[read_i:read_i + cigar_length],
                        length=cigar_length,
                        id=findZs("I"))]

    def mdDeletion():
        nonlocal md_i
        assert md[md_i] == '^'
        md_i += 1
        while md_i < len(md) and type(md[md_i]) is not int and md[md_i] in 'ACGT':
            md_i += 1
        return [Variant(typ="deletion",
                        ref=backbone,
                        pos=pos,
                        val=cigar_length,
                        length=cigar_length,
                        id=findZs("D"))]

    cmp_list = []
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


def readLink(index):
    """ Get related allele for the variant """
    id_alleles = {}
    for line in open(index + ".link"):
        id, alleles = line.strip().split('\t')
        id_alleles[id] = alleles.split()
    return id_alleles


def readVariants(index):
    """ Read hisat2 index.snp """
    id_alleles = readLink(index)
    variants = []
    for line in open(index + ".snp"):
        id, var_type, ref, pos, val = line.strip().split('\t')
        v = Variant(typ=var_type,
                    ref=ref,
                    pos=int(pos),
                    id=id,
                    allele=id_alleles.get(id, []),
                    val=val if var_type != "deletion" else int(val))
        variants.append(v)
    return variants


def readLocus(index):
    """ Read hisat2 xx.locus """
    # locus = 1-base position
    # variant = 0-base position
    # KIR2DL1*BACKBONE        KIR2DL1*BACKBONE        0       14761   14761   269-303 1267-1303 2031-2313 3754-4054 5588-5882 9039-9090 13360-13462 13924-13977 14075-14252 +
    gene_exons = {}
    for line in open(index + ".locus"):
        gene, _, _, _, _, exons, _ = line.split('\t')
        exons = map(lambda i:  list(map(lambda j: int(j) - 1, i.split("-"))), exons.split(' '))
        gene_exons[gene] = list(exons)
    return gene_exons


def isInExon(exons, variant):
    """ Is the variant in exon """
    # TODO: binary search, but naive method is not slow
    for exon in exons:
        if exon[0] <= variant.pos < exon[1]:
            return True
        if (variant.typ == "deletion" and \
            variant.pos < exon[0] and variant.pos + variant.val >= exon[0]):
            return True
    return False


def readAlleleLength(index):
    return {
        seq.id: len(seq.seq) for seq in SeqIO.parse(f"{index}_sequences.fa", "fasta")
    }


def findVariantId(variant):
    """ Match variant to hisat2 .snp """
    # The variant is not located in .snp.index but in .snp
    if variant in variants:
        # print("Find it", variant)
        return variants[variant]
    # cannot find: set it novel
    elif variant.typ in ["single", "insertion", "deletion"]:
        # print("Cannot find", variant)
        variant.id = f"nv{Variant.novel_count}"
        Variant.novel_count += 1
        return variant
    else:
        # print("Other", variant)
        assert variant.typ == "match"
        return variant


def extandMatch(variant_list):
    """ Remove novel variant """
    # force to transfer single novel variant to match
    v_list = []
    for v in variant_list:
        if v.typ == "single" and v.id.startswith("nv"):
            v.typ = "match"
        v_list.append(v)
    variant_list = v_list

    # merge match if previous type is also match
    v_list = []
    for v in variant_list:
        if v.typ == "match" and len(v_list) and v_list[-1].typ == "match":
            v_list[-1].length += v.length
        else:
            v_list.append(v)
    return v_list


def findAlts():
    pass
    # = typing_common.identify_ambigious_diffs(ref_seq,


def getVariantsBoundary(variant_list):
    """
    Find the lower and upper bound of list of variant

    Return [left, right) (of .snp variants)
    """
    left = variant_list[0].pos
    rigt = variant_list[-1].pos + variant_list[-1].length
    ref = variant_list[0].ref
    return (
        bisect.bisect_left(variants_sorted_list, Variant(ref=ref, pos=left, typ="single")),
        bisect.bisect_left(variants_sorted_list, Variant(ref=ref, pos=rigt, typ="single"))
    )


def getAlleleFromVariantList(variant_list, exon_only=False):
    """ Varinat -> list of allele """
    # Find all variant in the region
    left, right = getVariantsBoundary(variant_list)
    assert left <= right

    # remove novel or match
    variant_list = [v for v in variant_list if v.id and not v.id.startswith("nv")]

    # positive variation
    positive_var = set(v.id for v in variant_list)
    if exon_only:
        positive_allele = [v.allele for v in variant_list if v.in_exon]
    else:
        positive_allele = [v.allele for v in variant_list]

    # negative
    negative_allele = []
    for v in variants_sorted_list[left:right]:
        if v.id in positive_var:
            continue
        if exon_only and not v.in_exon:
            continue
        negative_allele.append(v.allele)
        # TODO: maybe some deletion is before left
        # they use gene_var_maxrights to deal with

    return positive_allele, negative_allele


def hisat2getCandidateAllele(positive_allele, negative_allele):
    candidate = None
    for allele in positive_allele:
        if candidate is None:
            candidate = set(allele)
        else:
            candidate &= set(allele)
    if candidate is None:
        return []

    for allele in negative_allele:
        candidate -= set(allele)
    return list(candidate)


def hisat2CountAllelePerPair(candidates):
    # count the allele
    # print(candidates)
    count = defaultdict(int)
    for allele in candidates:
        count[allele] += 1
    if not count:
        return set()

    # get maximum allele
    max_count = max(count.values())
    max_allele = set(i[0] for i in filter(lambda i: i[1] == max_count, count.items()))
    return max_allele


def hisatEMnp(allele_per_read):
    """
    single_abundance in hisat-genotype

    Accelerated version of EM - SQUAREM iteration
    Varadhan, R. & Roland, C. Scand. J. Stat. 35, 335-353 (2008)
    Also, this algorithm is used in Sailfish
    http://www.nature.com/nbt/journal/v32/n5/full/nbt.2862.html
    """
    # init probility (all allele has even probility)
    allele_name = sorted(set(allele for alleles in allele_per_read for allele in alleles))
    allele_map = dict(zip(allele_name, range(len(allele_name))))
    if norm_by_length:
        allele_len = np.array([seq_len[i] for i in allele_name])
    else:
        allele_len = np.ones(len(allele_name))
    allele_select_one = []
    for alleles in allele_per_read:
        x = np.zeros(len(allele_map))
        for i in map(allele_map.get, alleles):
            x[i] = 1
        allele_select_one.append(x)
    allele_select_one = np.array(allele_select_one)

    def getNextProb(prob_per_allele):
        a = prob_per_allele * allele_select_one
        b = a.sum(axis=1)[:, None]
        a = np.divide(a, b,
                      out=np.zeros(a.shape),
                      where=b!=0)
        a /= allele_len
        a = a.sum(axis=0)
        return a / a.sum()

    # Run EM
    prob = getNextProb(np.ones(len(allele_map)))
    for iters in range(0, iter_max):
        # SQUAREM
        prob_next = getNextProb(prob)
        prob_next2 = getNextProb(prob_next)
        r = prob_next - prob
        v = prob_next2 - prob_next - r
        r_sum = (r ** 2).sum()
        v_sum = (v ** 2).sum()
        if v_sum > 0.0:
            g = -np.sqrt(r_sum / v_sum)
            prob_next3 = prob - r * g * 2 + v * g ** 2
            prob_next3 = np.maximum(prob_next3, 0)
            prob_next = getNextProb(prob_next3)

        # Recalculate diff between previous stage
        diff = np.abs(prob - prob_next).sum()
        if diff <= diff_threshold:
            break

        # next
        # print(f"Iter: {iters} Diff: {diff}")
        prob = prob_next

    return dict(zip(allele_name, prob))


def hisat2StatAlleleCount(allele_counts, file=sys.stdout):
    # all candidates from all pairs
    count = defaultdict(int)
    for allele_count in allele_counts:
        for allele in allele_count:
            count[allele] += 1
    first_n = 10

    # sort and print
    count = sorted(count.items(), key=lambda i: i[1], reverse=True)
    for i, (allele, c) in enumerate(count[:first_n]):
        print(f"{i+1} {allele} (count: {c})", file=file)

    allele_prob = hisatEMnp(allele_counts)
    allele_prob = sorted(allele_prob.items(), key=lambda i: i[1], reverse=True)

    for i, (allele, prob) in enumerate(allele_prob[:first_n]):
        print(f"{i+1} ranked {allele} (abundance: {prob:.2f})", file=file)

    return count, allele_prob


def recordToVariants(record):
    variant_list, soft_clip = record2Variant(record)
    variant_list = map(findVariantId, variant_list)
    variant_list = extandMatch(variant_list)  # remove novel -> extend
    return variant_list


def writeBam(records, filename, filename_out):
    proc = subprocess.Popen(["samtools", "view", "-H", filename],
                            universal_newlines=True,
                            stdout=subprocess.PIPE,
                            stderr=open("/dev/null", 'w'))
    with open(filename_out, "w") as f:
        f.writelines(proc.stdout)
        f.writelines(records)


# index
def setIndex(index):
    global variants, variants_sorted_list, seq_len
    exon_region = readLocus(index)
    variants = readVariants(index)
    for v in variants:
        v.in_exon = isInExon(exon_region[v.ref], v)
    variants = {v: v for v in variants}
    variants_sorted_list = sorted([i for i in variants])
    seq_len = readAlleleLength(index)


def main(bam_file):
    new_suffix = ".hisatgenotype"
    # if exon_only: new_suffix += ".exon"
    name_out = os.path.splitext(bam_file)[0] + new_suffix

    # read bam file
    pair_reads = readPair(bam_file)
    pair_reads = filter(filterPair, pair_reads)

    save_reads = defaultdict(list)
    hisat_gene_alleles = defaultdict(list)
    tot_pair = 0
    for left_record, right_record in pair_reads:
        tot_pair += 1
        # if tot_pair > 100:
        #     break

        # group read by gene name
        backbone = left_record.split("\t")[2]

        lv = recordToVariants(left_record)
        rv = recordToVariants(right_record)

        # (left, right) x (positive, negative)
        lp, ln = getAlleleFromVariantList(lv)
        rp, rn = getAlleleFromVariantList(rv)

        # save the allele information for the read
        save_reads[backbone].append({
            'lp': lp, 'ln': ln, 'rp': rp, 'rn': rn,
            'l_sam': left_record,
            'r_sam': right_record,
        })

        # hisat2 method
        hisat_gene_alleles[backbone].append(hisat2CountAllelePerPair(
            hisat2getCandidateAllele(lp, ln) +
            hisat2getCandidateAllele(rp, rn)
        ))

        # TODO:
        # * error correction
        # * typing_common.identify_ambigious_diffs(ref_seq,
        # print(i, j)

    print("Filterd pairs", tot_pair)

    print(f"Save allele per reads in {name_out}.json")
    json.dump(dict(save_reads.items()), open(f"{name_out}.json", "w"))

    print(f"Save filtered bam in {name_out}.sam")
    records = []
    for reads in save_reads.values():
        for i in reads:
            records.append(i['l_sam'])
            records.append(i['r_sam'])
    writeBam(records, bam_file, name_out + ".sam")

    # hisat2 method
    hisat_result = []
    # hisat_report_f = open(f"{name_out}.report", "w")
    hisat_report_f = sys.stdout
    for backbone, hisat_alleles in sorted(hisat_gene_alleles.items(), key=lambda i: i[0]):
        print(backbone)
        hisat_result.append(hisat2StatAlleleCount(hisat_alleles, file=hisat_report_f))

    print(f"Save hisat-genotype calling in {name_out}.report*")
    json.dump(hisat_result, open(f"{name_out}.report.json", "w"))


if __name__ == "__main__":
    setIndex("index/kir_2100_raw.mut01")
    # name = "data/linnil1_syn_wide.00.kir_2100_raw.mut01.nosingle.bam"
    name = "data/linnil1_syn_wide.00.kir_2100_raw.mut01.bam"
    main(name)
