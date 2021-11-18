import re
import os
import pickle
import subprocess

import hisatgenotype_typing_common as typing_common
from hisatgenotype_typing_core import (
    read_Gene_vars_genotype_genome,
    read_backbone_alleles,
    read_Gene_alleles_from_vars,
    get_exonic_vars,
    get_rep_alleles
)


def get_bam_proc(alignment_fname, ref_locus, bam_sorted=True):
    # get bam file samtools command line
    alignview_cmd = ["samtools", "view", alignment_fname]
    _, chr, left, right = ref_locus[:4]
    alignview_cmd += ["%s:%d-%d" % (chr, left + 1, right + 1)]

    # samtools
    bamview_proc = subprocess.Popen(alignview_cmd,
                                    universal_newlines=True,
                                    stdout=subprocess.PIPE,
                                    stderr=open("/dev/null", 'w'))
    if bam_sorted:  # sort by name
        return bamview_proc
    sort_read_cmd = ["sort", "-k", "1,1", "-s"]  # -s for stable sorting
    alignview_proc = subprocess.Popen(sort_read_cmd,
                                      universal_newlines=True,
                                      stdin=bamview_proc.stdout,
                                      stdout=subprocess.PIPE,
                                      stderr=open("/dev/null", 'w'))
    return alignview_proc


def get_mpileup_cmd(alignment_fname, ref_locus):
    alignview_cmd = ["samtools", "view", alignment_fname]
    _, chr, left, right = ref_locus[:4]
    alignview_cmd += ["%s:%d-%d" % (chr, left + 1, right + 1)]
    return alignview_cmd


def readBam(alignment_fname, gene_vars, gene_var_list, ref_locus, ref_seq, base_locus=0):
    alignview_proc = get_bam_proc(alignment_fname, ref_locus, bam_sorted=False)
    mpileup_cmd    = get_mpileup_cmd(alignment_fname, ref_locus)
    mpileup        = typing_common.get_mpileup(mpileup_cmd, ref_seq, base_locus,
                                               gene_vars, False)  # , allow_discordant)
    return readBamCore(alignview_proc.stdout, gene_vars, gene_var_list, ref_seq, mpileup, base_locus)


def readBamCore(alignview_proc, gene_vars, gene_var_list, ref_seq, mpileup, base_locus=0):
    # parameters
    num_editdist     = 5
    allow_discordant = False
    verbose          = 0
    simulation       = False  # legency
    error_correction = False
    base_fname       = "kir"

    # setup
    cigar_re         = re.compile(r'\d+\w')  # Cigar regular expression
    left_read_ids    = set()
    right_read_ids   = set()
    var_count        = {}
    novel_var_count  = 0

    def debug_print(*arg):
        if verbose >= 2:
            print(*arg)

    # read line in bamfile
    for line in alignview_proc:
        line = line.strip()
        cols = line.split()
        read_id, flag, chr, pos, mapQ, cigar_str = cols[:6]

        node_read_id = read_id
        if simulation:
            read_id = read_id.split('|')[0]
        read_seq  = cols[9]
        read_qual = cols[10]
        flag      = int(flag)
        pos       = int(pos)
        pos      -= (base_locus + 1)
        if pos < 0:
            debug_print("Read is out of bound")
            continue

        # Unalined? Insurance that nonmaped reads will not be processed
        if flag & 0x4 != 0:
            if simulation:
                debug_print("Unaligned")
                debug_print("\t", line)
            continue

        # Concordantly mapped?
        if flag & 0x2 != 0:
            concordant = True
        else:
            concordant = False

        NM, Zs, MD, NH = "", "", "", ""
        for i in range(11, len(cols)):
            col = cols[i]
            if col.startswith("Zs"):
                Zs = col[5:]
            elif col.startswith("MD"):
                MD = col[5:]
            elif col.startswith("NM"):
                NM = int(col[5:])
            elif col.startswith("NH"):
                NH = int(col[5:])

        if NM > num_editdist:
            debug_print("Read > editdist")
            continue

        # Only consider unique alignment
        # linnil1: This doesn't affect alot
        if NH > 1:
            debug_print("Read not unique aligned")
            continue

        # Concordantly aligned mate pairs
        if not allow_discordant and not concordant:
            debug_print("Read is not concordant")
            continue

        # Add reads to nodes and assign left, right, or discordant
        is_left_read = flag & 0x40 != 0
        if is_left_read:            # Left read?
            if read_id in left_read_ids:
                debug_print("ERR")
                continue
            left_read_ids.add(read_id)
            if not simulation:
                node_read_id += '|L'
        elif flag & 0x80 != 0:      # Right read?
            if read_id in right_read_ids:
                debug_print("ERR")
                continue
            right_read_ids.add(read_id)
            if not simulation:
                node_read_id += '|R'
        else:
            assert allow_discordant
            if read_id in single_read_ids:
                debug_print("ERR")
                continue
            single_read_ids.add(read_id)
            if not simulation:
                node_read_id += '|U'

        if Zs:
            Zs_str = Zs
            Zs     = Zs.split(',')

        assert MD != ""
        MD_str_pos = 0
        MD_len     = 0
        Zs_pos     = 0
        Zs_i       = 0
        for _i in range(len(Zs)):
            Zs[_i]    = Zs[_i].split('|')
            Zs[_i][0] = int(Zs[_i][0])
        if Zs_i < len(Zs):
            Zs_pos += Zs[Zs_i][0]
        read_pos  = 0
        left_pos  = pos
        right_pos = left_pos
        cigars    = cigar_re.findall(cigar_str)
        cigars    = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
        cmp_list  = []
        num_error_correction = 0
        likely_misalignment  = False

        # Extract variants w.r.t backbone from CIGAR string
        softclip = [0, 0]
        for i in range(len(cigars)):
            cigar_op, length = cigars[i]
            if cigar_op == 'M':
                first       = True
                MD_len_used = 0
                cmp_list_i = len(cmp_list)
                while True:
                    if not first or MD_len == 0:
                        if MD[MD_str_pos].isdigit():
                            num = int(MD[MD_str_pos])
                            MD_str_pos += 1
                            while MD_str_pos < len(MD):
                                if MD[MD_str_pos].isdigit():
                                    num = num * 10 + int(MD[MD_str_pos])
                                    MD_str_pos += 1
                                else:
                                    break
                            MD_len += num
                    # Insertion or full match followed
                    if MD_len >= length:
                        MD_len -= length
                        if length > MD_len_used:
                            cmp_list.append(["match",
                                             right_pos + MD_len_used,
                                             length - MD_len_used])
                        break
                    first       = False
                    read_base   = read_seq[read_pos + MD_len]
                    MD_ref_base = MD[MD_str_pos]
                    MD_str_pos += 1
                    assert MD_ref_base in "ACGT"
                    if MD_len > MD_len_used:
                        cmp_list.append(["match",
                                         right_pos + MD_len_used,
                                         MD_len - MD_len_used])

                    _var_id = "unknown"
                    if read_pos + MD_len == Zs_pos and Zs_i < len(Zs):
                        assert Zs[Zs_i][1] == 'S'
                        _var_id = Zs[Zs_i][2]
                        Zs_i   += 1
                        Zs_pos += 1
                        if Zs_i < len(Zs):
                            Zs_pos += Zs[Zs_i][0]
                    else:
                        # Search for a known (yet not indexed)
                        # variant or a novel variant
                        ref_pos = right_pos + MD_len
                        var_idx = typing_common.lower_bound(gene_var_list,
                                                            ref_pos)
                        while var_idx < len(gene_var_list):
                            var_pos, var_id = gene_var_list[var_idx]
                            if var_pos > ref_pos:
                                break
                            if var_pos == ref_pos:
                                var_type, _, var_data = gene_vars[var_id]
                                if var_type == "single" \
                                        and var_data == read_base:
                                    _var_id = var_id
                                    break
                            var_idx += 1

                    cmp_list.append(["mismatch",
                                     right_pos + MD_len,
                                     1,
                                     _var_id])

                    MD_len_used = MD_len + 1
                    MD_len += 1
                    # Full match
                    if MD_len == length:
                        MD_len = 0
                        break

                # Correction for sequencing errors and
                # update for cmp_list
                if error_correction:
                    assert cmp_list_i < len(cmp_list)
                    name_readID = "aHSQ1008:175:C0JVFACXX:5:1109:17665:21583|L"
                    new_cmp_list, read_seq, _num_error_correction \
                        = error_correct(ref_seq,
                                        read_seq,
                                        read_pos,
                                        mpileup,
                                        gene_vars,
                                        gene_var_list,
                                        cmp_list[cmp_list_i:],
                                        node_read_id == name_readID)
                    cmp_list = cmp_list[:cmp_list_i] + new_cmp_list
                    num_error_correction += _num_error_correction

            elif cigar_op == 'I':
                _var_id = "unknown"
                if read_pos == Zs_pos and Zs_i < len(Zs):
                    assert Zs[Zs_i][1] == 'I'
                    _var_id = Zs[Zs_i][2]
                    Zs_i += 1
                    if Zs_i < len(Zs):
                        Zs_pos += Zs[Zs_i][0]
                else:
                    # Search for a known (yet not indexed)
                    # variant or a novel variant
                    var_idx = typing_common.lower_bound(gene_var_list,
                                                        right_pos)
                    while var_idx < len(gene_var_list):
                        var_pos, var_id = gene_var_list[var_idx]
                        if var_pos > right_pos:
                            break
                        if var_pos == right_pos:
                            var_type, _, var_data = gene_vars[var_id]
                            if var_type == "insertion" \
                                    and len(var_data) == length:
                                _var_id = var_id
                                break
                        var_idx += 1
                cmp_list.append(["insertion",
                                 right_pos,
                                 length,
                                 _var_id])
                if 'N' in read_seq[read_pos:read_pos+length]:
                    likely_misalignment = True

            elif cigar_op == 'D':
                if MD[MD_str_pos] == '0':
                    MD_str_pos += 1
                assert MD[MD_str_pos] == '^'
                MD_str_pos += 1
                while MD_str_pos < len(MD):
                    if not MD[MD_str_pos] in "ACGT":
                        break
                    MD_str_pos += 1
                _var_id = "unknown"
                if read_pos == Zs_pos and \
                   Zs_i < len(Zs) and \
                   Zs[Zs_i][1] == 'D':
                    _var_id = Zs[Zs_i][2]
                    Zs_i += 1
                    if Zs_i < len(Zs):
                        Zs_pos += Zs[Zs_i][0]
                else:
                    # Search for a known (yet not indexed) variant
                    # or a novel variant
                    var_idx = typing_common.lower_bound(gene_var_list,
                                                        right_pos)
                    while var_idx < len(gene_var_list):
                        var_pos, var_id = gene_var_list[var_idx]
                        if var_pos > right_pos:
                            break
                        if var_pos == right_pos:
                            var_type, _, var_data = gene_vars[var_id]
                            if var_type == "deletion" \
                                    and int(var_data) == length:
                                _var_id = var_id
                                break
                        var_idx += 1

                cmp_list.append(["deletion",
                                 right_pos,
                                 length,
                                 _var_id])

                # Check if this deletion is artificial alignment
                if right_pos < len(mpileup):
                    del_count, nt_count = 0, 0
                    for nt, value in mpileup[right_pos][1].items():
                        count = value[0]
                        if nt == 'D':
                            del_count += count
                        else:
                            nt_count += count

                    # DK - debugging purposes
                    if base_fname == "hla":
                        if del_count * 6 < nt_count:
                            likely_misalignment = True

            elif cigar_op == 'S':
                if i == 0:
                    softclip[0] = length
                    Zs_pos += length
                else:
                    assert i + 1 == len(cigars)
                    softclip[1] = length
            else:
                assert cigar_op == 'N'
                assert False
                cmp_list.append(["intron", right_pos, length])

            if cigar_op in "MND":
                right_pos += length

            if cigar_op in "MIS":
                read_pos += length

        # Remove softclip in cigar and modify read_seq and
        # read_qual accordingly
        if sum(softclip) > 0:
            if softclip[0] > 0:
                cigars = cigars[1:]
                read_seq = read_seq[softclip[0]:]
                read_qual = read_qual[softclip[0]:]
            if softclip[1] > 0:
                cigars = cigars[:-1]
                read_seq = read_seq[:-softclip[1]]
                read_qual = read_qual[:-softclip[1]]

            cigar_str = ""
            for type, length in cigars:
                cigar_str += str(length)
                cigar_str += type

        # if sum(softclip) > 0: #TODO Examine the purpose of this skip
        #     continue

        if right_pos > len(ref_seq):
            debug_print("Right pos is out of bound")
            continue

        if num_error_correction > max(1, num_editdist):
            debug_print("num_error_correction is larger")
            continue

        if likely_misalignment:
            debug_print("Read is misalign")
            continue

        # Add novel variants
        read_pos = 0
        for cmp_i in range(len(cmp_list)):
            type_, pos_, length_ = cmp_list[cmp_i][:3]
            if type_ != "match":
                var_id_ = cmp_list[cmp_i][3]
                if var_id_ == "unknown":
                    add = True
                    if type_ == "mismatch":
                        data_ = read_seq[read_pos]
                        if data_ == 'N':
                            add = False
                    elif type_ == "deletion":
                        data_ = str(length_)
                    else:
                        assert type_ == "insertion"
                        data_ = read_seq[read_pos:read_pos + length_]
                    if add:
                        if type_ != "mismatch":
                            type_add = type_
                        else:
                            type_add = "single"

                        var_id_, novel_var_count \
                            = add_novel_var(gene_vars,
                                            gene_var_list,
                                            novel_var_count,
                                            type_add,
                                            pos_,
                                            data_)
                        cmp_list[cmp_i][3] = var_id_
                if var_id_ != "unknown":
                    if var_id_ not in var_count:
                        var_count[var_id_] = 1
                    else:
                        var_count[var_id_] += 1

            if type_ != "deletion":
                read_pos += length_

        # cmp_list = [['match', 0, 171], ['mismatch', 171, 1, 'hv22'],
        #             ['match', 172, 90]]
        yield line, cmp_list
    # print(var_count, novel_var_count)
    # {'hv22': 62, 'nv0': 1, 'nv1': 1, 'nv2': 1, 'nv3': 1, 'nv4': 1, 'nv5': 1,
    #  'nv6': 1, 'hv37': 1, 'nv7': 1, 'nv8': 1, 'nv9': 1} 10


def add_novel_var(gene_vars,
                  gene_var_list,
                  novel_var_count,
                  var_type,
                  var_pos,
                  var_data):
    var_idx = typing_common.lower_bound(gene_var_list, var_pos)
    while var_idx < len(gene_var_list):
        pos_, id_ = gene_var_list[var_idx]
        if pos_ > var_pos:
            break
        if pos_ == var_pos:
            type_, _, data_ = gene_vars[id_]
            assert type_ != var_type or data_ != var_data
            if type_ != var_type:
                if var_type == "insertion":
                    break
                elif var_type == "single" and type_ == "deletion":
                    break
            else:
                if var_data < data_:
                    break
        var_idx += 1
    var_id = "nv%d" % novel_var_count
    assert var_id not in gene_vars
    gene_vars[var_id] = [var_type, var_pos, var_data]
    gene_var_list.insert(var_idx, [var_pos, var_id])
    return var_id, novel_var_count + 1


def read_hisat2(full_gg_path, base_fname):
    # Read locus
    refGenes       = {}
    refGene_loci   = {}
    typing_common.read_locus("%s.locus" % full_gg_path,
                             False,  # not genotype genome
                             base_fname,
                             refGenes,
                             refGene_loci)

    # Read .snp
    Vars, Var_list = read_Gene_vars_genotype_genome("%s.snp" % full_gg_path,
                                                    refGene_loci)

    # read .link
    Links          = typing_common.read_links("%s.link" % full_gg_path)
    # = {..., 'hv10087': ['KIR2DL1*001', 'KIR2DL1*002',...]}

    # Read .backbone
    Genes          = {}
    read_backbone_alleles(full_gg_path, refGene_loci, Genes)
    read_Gene_alleles_from_vars(Vars, Var_list, Links, Genes)

    return refGenes, refGene_loci, Vars, Var_list, Links, Genes


def get_allele_vars(gene_var_list, Links):
    allele_vars = {}
    for _, var_id in gene_var_list:
        if var_id not in Links:
            continue
        allele_list = Links[var_id]
        for allele_id in allele_list:
            # if allele_id not in Genes[gene]:
            #     continue
            if allele_id not in allele_vars:
                allele_vars[allele_id] = [var_id]
            else:
                allele_vars[allele_id].append(var_id)

    return allele_vars


def get_maxrights(gene_var_list, gene_vars):
    gene_var_maxrights = {}
    cur_maxright = -1
    for var_pos, var_id in gene_var_list:
        var_type, var_pos, var_data = gene_vars[var_id]
        if var_type == "deletion":
            var_pos = var_pos + int(var_data) - 1
        cur_maxright = max(cur_maxright, var_pos)
        gene_var_maxrights[var_id] = cur_maxright
    return gene_var_maxrights


def get_alts(ref_seq, allele_vars, gene_vars, gene_var_list):
    if False:  #os.path.exists("kir_alts.obj"):
        # read from python objects
        data = pickle.load(open("kir_alts.obj", "rb"))
        Alts_left = data['Alts_left']
        Alts_right = data['Alts_right']
        Alts_left_list = data['Alts_left_list']
        Alts_right_list = data['Alts_right_list']
    else:
        # calculate alts(slow)
        Alts_left, Alts_right = typing_common.get_alternatives(ref_seq,
                                                               allele_vars,
                                                               gene_vars,
                                                               gene_var_list,
                                                               False)
        Alts_left_list  = haplotype_alts_list(Alts_left, True)
        Alts_right_list = haplotype_alts_list(Alts_right, False)

        pickle.dump({'Alts_left': Alts_left,
                     'Alts_right': Alts_right,
                     'Alts_left_list': Alts_left_list,
                     'Alts_right_list': Alts_right_list},
                    open("kir_alts.obj", "wb"))

    return Alts_left, Alts_right, Alts_left_list, Alts_right_list


# For checking alternative alignments near the ends of alignments
def haplotype_alts_list(haplotype_alts, left=True):
    haplotype_list = []
    for haplotype in haplotype_alts.keys():
        if left:
            pos = int(haplotype.split('-')[-1])
        else:
            pos = int(haplotype.split('-')[0])
        haplotype_list.append([pos, haplotype])
    return sorted(haplotype_list, key=lambda x: x[0])


"""
# other init things

# alleles names
Gene_names = {}
for Gene_gene, data in Genes.items():
    Gene_names[Gene_gene] = list(data.keys())
# = {'KIR': ['KIR', "KIR3DS1*078"]}

# allele lengths
Gene_lengths = {}
for Gene_gene, Gene_alleles in Genes.items():
    Gene_lengths[Gene_gene] = {}
    for allele_name, seq in Gene_alleles.items():
        Gene_lengths[Gene_gene][allele_name] = len(seq)
# {'KIR': {'KIR': 21343}}
# Store de bruijn nodes that represent alleles
# allele_nodes = {}
# true_allele_nodes = {}

# Assembly graph
asm_graph = assembly_graph.Graph(ref_seq,
                                 gene_vars,
                                 ref_exons,
                                 ref_primary_exons,
                                 partial_alleles,
                                 true_allele_nodes,
                                 predicted_allele_nodes
                                 display_allele_nodes,
                                 simulation)

# Choose allele representives from those
# that share the primary exonic sequences
primary_exon_allele_reps, primary_exon_allele_rep_groups \
    = get_rep_alleles(Links, primary_exon_vars, allele_rep_set)
primary_exon_allele_rep_set \
    = set(primary_exon_allele_reps.values())

read_nodes    = []
read_var_list = []
"""
