#!/usr/bin/env python
"""
This file is modified from HISAT-genotype.
I rearrange and simplify it for HLA data. (2021/2/5 linnil1)

Feature:
* Find where the locus gene are in hg38
* Read IMGT
* Calculate consensus
* Write files for hisat2
* Write MSA to alignment format
* Build hisat2
"""

# --------------------------------------------------------------------------- #
# Copyright 2015, Daehwan Kim <infphilo@gmail.com>                            #
#                                                                             #
# HISAT-genotype is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# HISAT-genotype is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with HISAT-genotype.  If not, see <http://www.gnu.org/licenses/>.     #
# --------------------------------------------------------------------------- #

import os
from glob import glob
import sys
import subprocess
import re
import concurrent.futures
import copy

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from hisatgenotype_typing_common import (
    sort_genall, collapse_alleles
)
from hisatgenotype_typing_process import (
    split_haplotypes, key_varKey, hapKey)

from pyHLAMSA import HLAmsa, Genemsa

import pdb


# Main setting
base_fname         = "kir"          # This gene called hla
base_fullpath_name = "kir_merge"    # output name
consensus          = True           # Use consensus as backbone
threads            = 30             # Number of threads
collapse           = False          # Remove redundent allele. Turn off this when debugging

# basic setting
reference_index    = "./hg38.graph" # The hisat index built only with hg38.fa
fasta_dname        = "./IMGT/%s/fasta" % (base_fname.upper())
locus_list         = []             # gene name. e.g. A, B, DPA1
# inter_gap          = 30             # haplotype setting
# intra_gap          = 50             # haplotype setting
whole_haplotype    = True           # False -> Split allele
# partial            = True           # Merge nuc into gen
min_var_freq       = 0.01           # Filter snp with minimal freq into xx.index.snp

# some const
spliced_gene       = ['hla', 'rbg']
unspliced_gene     = ['codis', 'cyp', 'rrna']


def findBestAllele(gene):
    cigar_re = re.compile(r'\d+\w')
    aligner_cmd = ["hisat2"]
    # only include genes with no exact match to reference genome
    if base_fname not in ["cyp", "rbg"]:
        aligner_cmd += ["--score-min", "C,-12"]
    aligner_cmd += ["--no-unal",
                    "-x", reference_index,
                    "-f", "%s/%s_gen.fasta" % (fasta_dname, gene)]
    # print(aligner_cmd)
    align_proc = subprocess.Popen(aligner_cmd,
                                  universal_newlines=True,
                                  stdout=subprocess.PIPE,
                                  stderr=open("/dev/null", 'w'))
    allele_id   = ""
    best_chr    = ""
    best_left   = -1
    best_right  = -1
    best_AS     = -sys.maxsize
    best_strand = ""
    for line in align_proc.stdout:
        if line.startswith('@'):
            continue
        line = line.strip()
        cols = line.split()

        temp_allele_id, flag, chr, left, _, cigar_str = cols[:6]
        left   = int(left) - 1
        right  = left
        cigars = cigar_re.findall(cigar_str)
        cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
        if len(cigars) > 1 or cigars[0][0] != 'M':
            continue
        for i in range(len(cigars)):
            cigar_op, length = cigars[i]
            if cigar_op in "MND":
                right += length

        flag   = int(flag)
        strand = '-' if flag & 0x10 else '+'
        AS     = ""
        for i in range(11, len(cols)):
            col = cols[i]
            if col.startswith("AS"):
                AS = col[5:]
        assert AS != ""
        AS = int(AS)
        if AS > best_AS:
            allele_id   = temp_allele_id
            best_chr    = chr
            best_left   = left
            best_right  = right
            best_AS     = AS
            best_strand = strand

    align_proc.communicate()
    if allele_id == "":
        return best_chr, best_left, best_right, best_strand, "", ""

    # find allele name
    if base_fname in spliced_gene:
        allele_name = ""
        for line in open("%s/%s_gen.fasta" % (fasta_dname, gene)):
            line = line.strip()
            if not line.startswith('>'):
                continue
            if base_fname == 'hla':
                tmp_allele_id, tmp_allele_name = line[1:].split()[:2]
            else:
                tmp_allele_id   = line[1:].split()[0]
                tmp_allele_name = line[1:].split()[0]
            if allele_id == tmp_allele_id:
                allele_name = tmp_allele_name
                break
    else:
        allele_name = allele_id
    assert allele_name != "" and strand != ''

    print("%s-%s's reference allele is %s "
          "on '%s' strand of chromosome %s from %d %d" % (
              base_fname.upper(), gene, allele_name,
              strand, chr, best_left, best_right),
          file=sys.stderr)
    assert best_chr != "" and best_left >= 0 and best_right > best_left
    return best_chr, best_left, best_right, best_strand, allele_id, allele_name


def count2freqdict(consensus_freq):
    # This function is part of create_consensus_seq
    # Convert a list form of consensus_freq to a dictionary form
    temp_freq = []
    for j in range(len(consensus_freq)):
        freq_dic = {}
        s = sum(consensus_freq[j])
        for k in range(len(consensus_freq[j])):
            freq = consensus_freq[j][k]
            if freq == 0:
                continue
            nt = "ATCG."[k]
            freq_dic[nt] = freq / s
        temp_freq.append(freq_dic)
    consensus_freq = temp_freq
    return consensus_freq


def checkOmit(backbone_seq):
    # Check for empty sequences or omitted nucleotides.
    # This can be empty if there is a problem with the database.
    # Solution is to omit the offending gene from the
    # database to prevent premature termination of database build.
    omits = ['.', 'E', '~']
    breakout = False

    for omit in omits:
        if omit in backbone_seq:
            if verbose:
                print('%s in backbone of %s with no "\
                          "minimum variation set' % (omit, gene))
            print('Error in database: Omitting %s!!' % (gene),
                  file=sys.stderr)
            breakout = True
            break

    if breakout:
        raise "Omit"


def msa2Var(backbone_seq, backbone_freq, names, seqs):
    # Insert Variants into Var dictionary variable
    # Inherits 4 variables from above and shouldn't be moved
    def insertVar(type, info):
        pos, backbone_pos, data = info
        if type in "MI":
            varKey = "%d-%s-%s" % (pos, type, data)
        else:
            varKey = "%d-%s-%d" % (pos, type, data)

        if varKey not in Vars:
            if type == 'M':
                assert data in backbone_freq[backbone_pos], \
                    "Data %s not in backbone %s of type %s" % (
                        data,
                        backbone_freq[backbone_pos],
                        type)
                freq = backbone_freq[backbone_pos][data]
            elif type == 'D':
                del_len = int(data)
                freq = 100.0
                for d in range(del_len):
                    assert '.' in backbone_freq[backbone_pos + d]
                    freq2 = backbone_freq[backbone_pos + d]['.']
                    if freq2 < freq:
                        freq = freq2
            else:
                assert type == 'I'
                ins_len = len(data)
                freq = 100.0
                acc = 0
                for i in range(ins_len):
                    nt = data[i]
                    # TOFIX
                    # both query and backbone are '.'
                    while nt not in backbone_freq[backbone_pos + i + acc] \
                            and backbone_freq + i + acc < len(backbone_freq):
                        acc += 1
                    assert backbone_pos + i + acc < len(backbone_freq)
                    # if nt not in backbone_freq[backbone_pos + i + acc]:
                    #    print("NT not exist", info, i, nt, backbone_freq[backbone_pos + i + acc])
                    # assert nt in backbone_freq[backbone_pos + i + acc]
                    freq2 = backbone_freq[backbone_pos + i + acc][nt]
                    if freq2 < freq:
                        freq = freq2
                if not consensus:
                    assert freq <= min_var_freq
            Vars[varKey] = [freq, [cmp_name]]
        else:
            Vars[varKey][1].append(cmp_name)

    # nt -> var
    Vars = {}
    seq_len = len(backbone_seq)
    for cmp_name, id in names.items():
        # if cmp_name == backbone_name:
        #     continue
        assert id < len(seqs)
        cmp_seq = seqs[id]
        if len(cmp_seq) != seq_len:
            print("Warning: the length of %s (%d) is different from %d" % (
                  cmp_name, len(cmp_seq), seq_len),
                  file=sys.stderr)
            continue

        insertion = []
        deletion = []
        ndots = 0
        for s in range(seq_len):
            assert not (insertion and deletion)
            bc = backbone_seq[s]
            cc = cmp_seq[s]
            # print(s, bc, cc, backbone_freq[s])  # debug when assert error
            assert cc in backbone_freq[s]
            if bc not in '.~' and cc not in '.~':
                if insertion:
                    insertVar('I', insertion)
                    insertion = []
                elif deletion:
                    insertVar('D', deletion)
                    deletion = []
                if bc != cc:
                    mismatch = [s - ndots, s, cc]
                    insertVar('M', mismatch)
            elif bc == '.' and cc not in '.~':
                if deletion:
                    insertVar('D', deletion)
                    deletion = []
                if insertion:
                    insertion[2] += cc
                else:
                    insertion = [s - ndots, s, cc]
            elif bc not in '.~' and cc == '.':
                if insertion:
                    insertVar('I', insertion)
                    insertion = []
                if deletion:
                    deletion[2] += 1
                else:
                    deletion = [s - ndots, s, 1]

            if bc == '.':
                ndots += 1

        if insertion:
            insertVar('I', insertion)
        elif deletion:
            insertVar('D', deletion)
    return Vars


def var2Allele(Vars):
    # var_name to key
    Vars_ = {}
    for key, values in Vars.items():
        freq, names_ = values
        for name in names_:
            if name not in Vars_:
                Vars_[name] = [key]
            else:
                Vars_[name].append(key)
    for name, vars in Vars_.items():
        Vars_[name] = sorted(vars, key=key_varKey)

    return Vars_


def getExonStr(backbone_seq, msa):
    exon_count = []
    pos = msa._calculate_position()

    for i in range(1, len(msa.blocks), 2):
        exon_count.append(
            len(msa.select_chunk([i]).select_complete().alleles))

    exon_str = ""
    for i in range(1, len(msa.blocks), 2):
        exon_left  = len(backbone_seq[:pos[i + 0]].replace(".", "").replace("E", ""))
        exon_right = len(backbone_seq[:pos[i + 1]].replace(".", "").replace("E", ""))
        primary = exon_count[i // 2] == max(exon_count)
        if exon_str != "":
            exon_str += ","
        exon_str += "%d-%d%s" % (exon_left + 1, exon_right + 1, 'p' if primary else '')
    return exon_str


def groupHaplo(Vars):
    Vars_ = var2Allele(Vars)

    # exclude snp < min_freq
    keys = sorted(Vars.keys(), key=key_varKey)
    excluded_vars = set()
    for key in keys:
        locus, type, data = key.split('-')
        if type == "M" and Vars[key][0] < min_var_freq:
            excluded_vars.add(key)

    # group haplotype
    all_haplotypes = []
    i = 0
    while i < len(keys):
        var_leftmost = sys.maxsize
        var_rightmost = -1

        # the haplo is from i to j
        key_i = keys[i]
        locus, type, data = key_i.split('-')
        locus = int(locus)
        if type == 'D':
            locus += (int(data) - 1)
        prev_locus = locus
        if whole_haplotype:
            j = len(keys)
        else:
            j = i + 1
            while j < len(keys):
                key_j = keys[j]
                locus2, type2, data2 = key_j.split('-')
                locus2 = int(locus2)
                if prev_locus + inter_gap < locus2:
                    break
                prev_locus = locus2
                if type == 'D':
                    prev_locus += (int(data) - 1)
                j += 1

        # get available allele in i to j
        alleles = set()
        for key_k in keys[i:j]:
            freq, names_ = Vars[key_k]
            locus, type, data = key_k.split('-')
            if type == "M" and freq < min_var_freq:
                continue
            alleles |= set(names_)

            # Update leftmost and rightmost of Vars
            locus, type, data = key_k.split('-')
            left = right = int(locus)
            if type == 'D':
                right = left + int(data) - 1
            if left < var_leftmost:
                var_leftmost = left
            if var_rightmost < right:
                var_rightmost = right
        assert var_leftmost <= var_rightmost

        # group haplotype for each allele in i to j
        haplotypes = set()
        cur_vars = set(keys[i:j]) - excluded_vars
        for allele in alleles:
            allele_vars = set(Vars_[allele]) - excluded_vars
            allele_var_sort = sorted(list(cur_vars & allele_vars),
                                     key=key_varKey)

            allele_cur_vars = '#'.join(allele_var_sort)
            haplotypes.add(allele_cur_vars)

        if not whole_haplotype:
            haplotypes = split_haplotypes(haplotypes, intra_gap)

        haplotypes = sorted(list(haplotypes), key=hapKey)
        all_haplotypes.append([haplotypes, var_leftmost, var_rightmost])

        # next
        i = j

    return all_haplotypes


def gene2hisat(gene, msa):
    msa = msa.shrink()
    print(msa.blocks)
    print(msa.labels)

    print(f"[{gene}] Calculate consensus", file=sys.stderr)
    backbone_freq = msa.calculate_frequency()
    if f"{gene}*BACKBONE" not in msa.alleles:
        msa.add(f"{gene}*BACKBONE", msa.get_consensus(include_gap=False))

    backbone_seq = msa.get(f"{gene}*BACKBONE")
    backbone_name = f"{gene}*BACKBONE"

    backbone_freq = count2freqdict(backbone_freq)
    msa_ori = copy.deepcopy(msa)        # for writing msa
    msa.fill_imcomplete(backbone_name)  # for building hisat2
    print(f"[{gene}] Backbone len={len(backbone_seq)} and " + \
          f"no-gap len={len(backbone_seq.replace('-', ''))}")

    print(f"[{gene}] Find and remove collapse seqs",
          file=sys.stderr)
    # msa -> hisatmsa,  "-" -> "."
    names = {}
    seqs = []
    for name, seq in msa.alleles.items():
        if name != backbone_name:
            names[name] = len(seqs)
            seqs.append(seq.replace("-", "."))

    if collapse:
        names, seqs, collapsed = collapse_alleles(names,
                                                  seqs,
                                                  list_collapse=True,
                                                  verbose=True)
    # hisatmsa -> msa
    msa.alleles = {}
    msa.add(backbone_name, backbone_seq)  # no changed
    for name in names:
        msa.alleles[name] = seqs[names[name]].replace(".", "-")

    # msa -> hisatmsa
    backbone_seq = backbone_seq.replace("-", ".")

    # Code below mostly same in hisat-genotype
    # check backbone did not contain gap when min_freq = 0
    if min_var_freq <= 0.0:
        checkOmit(backbone_seq)

    print(f"[{gene}] Variant to key string", file=sys.stderr)
    Vars = msa2Var(backbone_seq, backbone_freq, names, seqs)
    print(f"[{gene}] Number of variants is {len(Vars)}.",
          file=sys.stderr)

    print(f"[{gene}] Find haplotype", file=sys.stderr)
    haplotypes = groupHaplo(Vars)

    # format locus info
    exon_str = getExonStr(backbone_seq, msa_ori)
    out_locus = ""

    allele_id = None
    # chr, left, right, strand, allele_id, allele_name = findBestAllele(gene)
    if allele_id is None:
        chr, left, right, strand = f"{gene}*BACKBONE", 0, len(backbone_seq) - 1, "+"
    out_locus = "%s\t%s\t%d\t%d\t%d\t%s\t%s" % (
        backbone_name, chr,
        left, right - 1,
        len(backbone_seq.replace('.', '')),
        exon_str,
        strand)

    print(f"[{gene}] Done", file=sys.stderr)
    return ([backbone_name, msa, msa_ori],
            out_locus,
            [backbone_name, Vars, haplotypes])


def writeSeq(backbone_file, name, seq):
    assert "~" not in seq
    SeqIO.write(SeqRecord(Seq(seq), id=name, description=""),
                backbone_file,
                "fasta")


def write(out_msa, out_locus, out_vars):
    # Open file
    locus_file     = open(base_fullpath_name + ".locus", 'w')
    # Write the backbone sequences into a fasta file
    backbone_file  = open(base_fullpath_name + "_backbone.fa", 'w')
    # variants w.r.t the backbone sequences into a SNP file
    var_file       = open(base_fullpath_name + ".snp", 'w')
    var_index_file = open(base_fullpath_name + ".index.snp", 'w')
    # variant frequence
    var_freq_file  = open(base_fullpath_name + ".snp.freq", 'w')
    # haplotypes
    haplotype_file = open(base_fullpath_name + ".haplotype", 'w')
    # pairs of a variant and the corresponding HLA alleles into a LINK file
    link_file      = open(base_fullpath_name + ".link", 'w')
    # Write all the sequences with dots removed into a file
    input_file     = open(base_fullpath_name + "_sequences.fa", 'w')
    # Write allele names into a file
    allele_file    = open(base_fullpath_name + ".allele", 'w')
    # Read partial alleles from hla.data, and write them into a file
    partial_file   = open(base_fullpath_name + ".partial", 'w')

    print(f"Write locus", file=sys.stderr)
    for locus_str in out_locus:
        print(locus_str, file=locus_file)

    # Write
    #       (1) variants w.r.t the backbone sequences into a SNP file
    #       (2) pairs of a variant and the corresponding HLA alleles
    #               into a LINK file
    num_vars = 0
    num_haplotypes = 0
    for backbone_name, Vars, all_haplotypes in out_vars:
        print(f"Write {backbone_name} snp(link snp index.snp freq)",
              file=sys.stderr)
        keys = sorted(Vars.keys(), key=key_varKey)
        var2ID = {}
        for key in keys:
            freq, names_ = Vars[key]
            names_ = sorted(names_)
            varID = "hv%d" % (num_vars)

            # format
            locus, type, data = key.split('-')
            locus = int(locus)
            if type == 'M':
                type_str = "single"
            elif type == 'I':
                type_str = "insertion"
            else:
                assert type == 'D'
                type_str = "deletion"
            base_locus = 0
            tmp_backbone_name = backbone_name
            var_str = "%s\t%s\t%s\t%d\t%s" \
                % (varID,
                   type_str,
                   tmp_backbone_name,
                   base_locus + locus,
                   data)

            # write to file
            print(var_str, file=var_file)
            if type != "M" or freq >= min_var_freq:
                print(var_str, file=var_index_file)
            print("%s\t%.2f" % (varID, freq),
                  file=var_freq_file)
            print("%s\t%s" % (varID, ' '.join(names_)),
                  file=link_file)
            var2ID[key] = num_vars
            num_vars += 1

        # Write haplotypes
        print(f"Write {backbone_name} haplotypes", file=sys.stderr)
        add_seq_len = 0
        for haplotypes, var_leftmost, var_rightmost in all_haplotypes:
            assert var_leftmost <= var_rightmost
            # sanity_vars = set()
            for h_i in range(len(haplotypes)):
                h = haplotypes[h_i].split('#')
                varIDs = []
                for var in h:
                    varIDs.append("hv%s" % var2ID[var])
                    # sanity_vars.add(var2ID[var])
                h_begin = var_leftmost
                h_end   = var_rightmost
                base_locus = 0
                tmp_backbone_name = backbone_name
                print("ht%d\t%s\t%d\t%d\t%s" % (
                    num_haplotypes,
                    tmp_backbone_name,
                    base_locus + h_begin,
                    base_locus + h_end,
                    ','.join(varIDs)),
                    file=haplotype_file)
                num_haplotypes += 1
                add_seq_len += (h_end - h_begin + 1)
            # assert len(sanity_vars) == len(cur_vars)

        print("Length of additional sequences for haplotypes:", add_seq_len,
              file=sys.stderr)

    # Write all the sequences with dots removed into a file
    print(f"Write backbone and partial and allele and input", file=sys.stderr)
    for backbone_name, msa, msa_ori in out_msa:
        # Backbone
        writeSeq(backbone_file, backbone_name,
                 msa.get(backbone_name).replace('-', ''))

        # sorted_all_alleles = sort_genall(list(msa.alleles.keys()),
        #                                  alleles=True)
        partial_alleles = msa_ori.select_complete().alleles

        for name in msa.alleles.keys():
            # Seq
            if name == backbone_name:
                continue
            writeSeq(input_file, name,
                     msa.get(name).replace("E", "").replace("-", ""))
            # Allele name
            print(name, file=allele_file)

            # exon-only allele
            if name not in partial_alleles:
                print(name, file=partial_file)

        open(f"{base_fullpath_name}.{msa.gene_name}_gen.txt", 'w').write(
            msa_ori.select_complete().format_alignment_diff(backbone_name))
        open(f"{base_fullpath_name}.{msa.gene_name}_nuc.txt", 'w').write(
            msa_ori.select_exon().format_alignment_diff(backbone_name))
        # already reverse the label when - strand
        msa.save_gff(f"{base_fullpath_name}.{msa.gene_name}.gff", strand="+")
        msa.save_bam(f"{base_fullpath_name}.{msa.gene_name}.bam", backbone_name)

    # close file
    print(f"Done", file=sys.stderr)
    backbone_file.close()
    locus_file.close()
    var_file.close()
    var_index_file.close()
    var_freq_file.close()
    haplotype_file.close()
    link_file.close()
    input_file.close()
    allele_file.close()
    partial_file.close()


def main(name):
    # default
    global base_fullpath_name
    base_fullpath_name = name

    pool = concurrent.futures.ProcessPoolExecutor(threads)

    # linnil1: another variable?
    if base_fullpath_name.startswith("kir_merge"):
        msa = Genemsa.load_msa(f"{base_fullpath_name}.save.fa", f"{base_fullpath_name}.save.gff")
        msa.seq_type = "gen"
        msa.gene_name = base_fullpath_name
        pool_gene = [gene2hisat("KIR", msa)]
    elif base_fullpath_name.startswith("kir_split"):
        merge_name = base_fullpath_name.replace("split", "merge")
        msa = Genemsa.load_msa(f"{merge_name}.save.fa", f"{merge_name}.save.gff")
        msa.seq_type = "gen"
        msa.gene_name = "kir_merge"
        pool_gene = []
        for gene in set(i.split("*")[0] for i in msa.alleles.keys() if "BACKBONE" not in i):
            newmsa = msa.select_allele(gene + ".*")
            newmsa.gene_name = gene
            pool_gene.append(gene2hisat(gene, newmsa))
        """
        pool_gene = []
        for i in glob("kir_split.*.save.fa"):
            gene = i.split('.')[1]
            msa = Genemsa.load_msa(i, i[:-3] + ".gff")
            msa.seq_type = "gen"
            msa.gene_name = gene
            pool_gene.append(gene2hisat(gene, msa))
        """
    else:
        print("Not kir_merge or kir_split", file=sys.stderr)
        sys.exit()

    # write file
    out_msa = []
    out_locus = []
    out_vars = []
    for gene in pool_gene:
        msa, locus, vars = gene
        if len(msa):
            out_msa.append(msa)
        if len(locus):
            out_locus.append(locus)
        if len(vars):
            out_vars.append(vars)

    write(out_msa, out_locus, out_vars)
    os.system(f"ln -s {base_fullpath_name}_backbone.fa {base_fullpath_name}.fa")


def build(name):
    global base_fullpath_name
    base_fullpath_name = name
    cmd = f"""\
        hisat2-build {base_fullpath_name}_backbone.fa \
                     --snp {base_fullpath_name}.index.snp \
                     --haplotype {base_fullpath_name}.haplotype \
                     -p {threads} \
                     --verbose {base_fullpath_name}.graph
    """
    print(cmd)
    os.system(cmd)


# dk -v $PWD/hisat-genotype:/root/hisatgenotype linnil1/hisat2
# apt update
# apt install python3-pip
# pip3 install biopython
# pip3 install -e pyHLAMSA
if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "kir_merge"
    main(name)
    build()
