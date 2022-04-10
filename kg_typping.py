# This script is copy from hisat-genotype/hisatgenotype_typing_core
import re
import subprocess

num_editdist = 5
NH_filtering = False


def readBam(alignment_fname):
    # get bam file samtools command line
    alignview_cmd = ["samtools", "view", "-h", alignment_fname,
                     "|", "samtools", "sort", "-n", "-", "-O", "SAM"]
    proc = subprocess.Popen(" ".join(alignview_cmd),
                            shell=True,
                            universal_newlines=True,
                            stdout=subprocess.PIPE,
                            stderr=open("/dev/null", 'w'))
    return proc.stdout


def readPair(alignment_fname):
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
            print("Skip", line)
            continue

        # If pair read in temp dict -> remote from temp, check and yield
        num_reads += 1
        next_ref = ref
        next_id = (next_ref, next_pos)
        if next_id in reads:
            next_line = reads[next_id]
            # is left and right
            if ((int(next_line.split('\t')[1]) | int(flag)) & (64 + 128)) != 64 + 128:
                print("Strange case", line, next_line)
                continue
            # main
            del reads[next_id]
            num_pairs += 1
            yield line, next_line
        # if not find -> save in temp
        else:
            reads[(ref, pos)] = line

    # summary
    print("Reads:", num_reads, "Pairs:", num_pairs)


def filterRead(line):
    flag = int(line.split('\t')[1])

    # Concordantly mapped?
    if flag & 2 == 0:
        print("Not concordantly mapped read")
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
        print("Read not unique aligned")
        return False
    if NM > num_editdist:
        print("Read > editdist")
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
            # remove 0
            md = filter(None, md)

    return list(zs), list(md)


def record2Variant(line):
    """
    Turn bam record into a list of variant format

    Variant format = (type, pos, length, varid, seq_in_that_position)
    * type: match, mismatch, insertion, deletion
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
                yield ["match", pos + md_len_last, cigar_length - md_len_last]
                break

            # read the base
            # * md[md_i] is reference
            # * read_seq[pos] is alt
            read_base = read_seq[read_i + md_len]
            assert md[md_i] in "ACGT"
            assert md[md_i] != read_base
            md_i += 1

            # the sequence before current_pos is the same
            if md_len > md_len_last:
                yield ["match", pos + md_len_last, md_len - md_len_last]

            # the current position is mismatch
            yield ["mismatch", pos + md_len, 1, findZs("S"), read_base]
            md_len += 1
            md_len_last = md_len

            # cigar OK, increase the pos and set md_len = md_len_last = 0
            if md_len == cigar_length:
                md_len = 0
                break

    def mdInsertion():
        # md does not contains insertion info
        return [["insertion", pos, cigar_length, findZs("I")]]

    def mdDeletion():
        nonlocal md_i
        assert md[md_i] == '^'
        md_i += 1
        while md_i < len(md) and type(md[md_i]) is not int and md[md_i] in 'ACGT':
            md_i += 1
        return [["deletion", pos, cigar_length, findZs("D")]]

    cmp_list = []
    soft_clip = [0, 0]
    # print(line)
    for cigar_i, (cigar_op, cigar_length) in enumerate(cigars):
        # print("cigar", cigar_op, cigar_length)
        # print("pos read_i", pos, read_i)
        # print("ZS", zs, zs_i, zs_pos)
        # print('MD', md, md_i, md_len)
        # print('cmp_list', cmp_list)
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

    # double check
    assert zs_i == len(zs)
    assert md_i == len(md)
    assert read_i == len(read_seq)

    return cmp_list, soft_clip


# out = readBam("data/linnil1_syn_wide.00.kir_2100_raw.mut01.nosingle.bam")
out = readPair("data/linnil1_syn_wide.00.kir_2100_raw.mut01.nosingle.bam")
tot_pair = 0
for i, j in filter(filterPair, out):
    tot_pair += 1
    record2Variant(i)
    record2Variant(j)
    # TODO:
    # * error correction
    # * find var_id in .snp
    # print(i, j)
print("Filterd pairs", tot_pair)
