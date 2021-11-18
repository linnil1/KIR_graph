from collections import defaultdict
from pprint import pprint
from copy import deepcopy
import os
import subprocess
import json
import glob
import sys

# hisat
# sys.path.insert(0, "/home/linnil1/hisat2/hisat-genotype/hisatgenotype_modules/")
import hisatgenotype_typing_common as typing_common
import hisatgenotype_assembly_graph as assembly_graph

# local
from pipeline_typing_util import (
    readBam, read_hisat2,
    get_allele_vars, get_maxrights, get_alts,
    get_exonic_vars, get_rep_alleles
)

base_fname          = "kir"  # not needed to change
allow_nonpair       = False


def hisatdataInit(path):
    global full_gg_path
    full_gg_path = path

    print("Reading hisat data")
    global refGenes, refGene_loci, Vars, Var_list, Links, Genes, Alts
    refGenes, refGene_loci, Vars, Var_list, Links, Genes \
                        = read_hisat2(full_gg_path, base_fname)

    print("Alts")
    Alts = {}
    """
    for gene, ref_allele in refGenes.items():
        Alts[gene] = [[],[],[],[]]
    """
    for gene, ref_allele in refGenes.items():
        Alts[gene] = get_alts(Genes[gene][ref_allele],
                              get_allele_vars(Var_list[gene], Links),
                              Vars[gene],
                              Var_list[gene])
    print("Alts Done")


def readPair(*args):
    num_reads      = 0
    num_pairs      = 0
    prev_read_id   = None
    prev_read_lr   = 0
    prev_line_l    = None
    prev_line_r    = None
    prev_cmp_l     = None
    prev_cmp_r     = None

    for line, cmp_list in readBam(*args):
        num_reads += 1
        read_id, flag = line.split()[:2]

        # Clear
        if read_id != prev_read_id:
            prev_read_id  = read_id
            prev_read_lr  = 0

        # save left or right read
        is_left_read = int(flag) & 0x40 != 0
        if is_left_read:
            prev_read_lr |= 1
            prev_line_l   = line
            prev_cmp_l    = cmp_list
        else:
            prev_read_lr |= 2
            prev_line_r   = line
            prev_cmp_r    = cmp_list

        if prev_read_lr == 3:
            num_pairs += 1
            yield prev_line_l, prev_line_r, prev_cmp_l, prev_cmp_r

    # summary
    print("Reads:", num_reads, "Pairs:", num_pairs)


def hisatTyping(gene):
    # -----------------------
    # --- Load Hisat data ---
    # -----------------------
    ref_allele          = refGenes[gene]
    # = "KIR"
    ref_locus           = refGene_loci[gene]
    # = ['KIR', 'KIR*BACKBONE', 0, 21342, [[551, 591], [3047, 3083], [3994, 4279], [5785, 6085], [7737, 8031], [11544, 11595], [19638, 19744], [20207, 20260], [20380, 20654]], [[551, 591], [3047, 3083], [3994, 4279], [5785, 6085], [7737, 8031], [11544, 11595], [19638, 19744], [20207, 20260], [20380, 20654]]]
    base_locus          = ref_locus[2]
    # ref_exons           = ref_locus[-2]
    # ref_primary_exons   = ref_locus[-1]

    gene_vars           = deepcopy(Vars[gene])
    # = {... 'hv10086': ['single', 21332, 'C'], 'hv10087': ['deletion', 21336, '7']}
    gene_var_list       = deepcopy(Var_list[gene])
    # = [... [21332, 'hv10086'], [21336, 'hv10087']]

    ref_seq             = Genes[gene][ref_allele]
    gene_names          = Genes[gene].keys()
    # ref_seq = ATCG...
    # ['KIR', 'KIR3DS1*078', 'KIR3DS1*013']

    # List of nodes that represent alleles
    allele_vars         = get_allele_vars(gene_var_list, Links)
    # = {... 'KIR2DL5A*034': ['hv23', 'hv54', 'hv86', ...]}

    # Extract variants that are within exons
    # exon_vars           = get_exonic_vars(gene_vars, ref_exons)
    # primary_exon_vars   = get_exonic_vars(gene_vars, ref_primary_exons)

    # Choose allele representives from those
    # that share the same exonic sequences
    # allele_reps, allele_rep_groups \
    #                     = get_rep_alleles(Links, exon_vars)
    # {... 'KIR2DS2*021': 'KIR2DS2*021'}
    # {... 'KIR2DS2*021': ['KIR2DS2*021']}
    # allele_rep_set      = set(allele_reps.values())
    # {'KIR3DL3*108', ...}

    gene_var_maxrights  = get_maxrights(gene_var_list, gene_vars)
    # {..., 'hv10087': 21342}

    Alts_left, Alts_right, Alts_left_list, Alts_right_list = Alts[gene]
    # {...21283-hv10066-hv10067-hv10069-21286: {'21283-hv10064-hv10067-hv10069-21287'}
    # [...[8950, '8821-hv4802-hv4803-hv4804-8950']]

    # --------------------------------
    # --- hisat-kir Functions --------
    # ----- Modified from add_count
    # --------------------------------
    def count_allele_pn(ht):
        """
        Transfer ht to allele by .link file and using AND

        params:
          ht: example: 5366-hv2610-hv2659-hv2683-5668, 19309-19609
        """
        orig_ht = ht
        ht = ht.split('-')
        assert len(ht) >= 2
        left  = int(ht[0])
        right = int(ht[-1])
        assert left <= right

        # use full name
        def getName(a):
            return a

        # read locus of each base shown in bamfile
        ht = ht[1:-1]
        for i in range(len(ht)):
            var_id = ht[i]
            # remove novel
            if var_id.startswith("nv") or var_id not in Links:
                continue
            # positive and negative allele
            for i in getName(Links[var_id]):
                allele_count_p[i] += 1
            for i in getName(set(gene_names) - set(Links[var_id])):
                allele_count_n[i] += 1
        ht = set(ht)

        # read locus of each base NOT shown in bamfile but in index
        tmp_alleles = set()
        var_idx = typing_common.lower_bound(gene_var_list, right + 1)
        var_idx = min(var_idx, len(gene_var_list) - 1)
        while var_idx >= 0:
            _, var_id = gene_var_list[var_idx]
            # remove novel
            if var_id.startswith("nv") \
                    or var_id in ht \
                    or var_id not in Links:
                var_idx -= 1
                continue
            if var_id in gene_var_maxrights \
                    and gene_var_maxrights[var_id] < left:
                break
            var_type, var_left, var_data = gene_vars[var_id]
            var_right = var_left
            if var_type == "deletion":
                var_right = var_left + int(var_data) - 1
            if (var_left >= left and var_left <= right) \
                    or (var_right >= left and var_right <= right):
                tmp_alleles |= set(Links[var_id])
                # positive and negative allele
                for i in getName(Links[var_id]):
                    allele_count_n[i] += 1
                for i in getName(set(gene_names) - set(Links[var_id])):
                    allele_count_p[i] += 1
            var_idx -= 1

    # ------------------------
    # --- Hisat2 Functions ---
    # ------------------------
    def add_count(count_per_read, ht, add):
        if base_fname == "genome" and len(count_per_read) == 1:
            for allele in count_per_read.keys():
                count_per_read[allele] = add
            return

        orig_ht = ht
        ht = ht.split('-')

        assert len(ht) >= 2
        left  = int(ht[0])
        right = int(ht[-1])
        assert left <= right

        ht = ht[1:-1]
        alleles = set(Genes[gene].keys()) - set([ref_allele])
        for i in range(len(ht)):
            var_id = ht[i]
            if var_id.startswith("nv") or \
               var_id not in Links:
                continue
            alleles &= set(Links[var_id])
        ht = set(ht)

        tmp_alleles = set()
        var_idx = typing_common.lower_bound(gene_var_list, right + 1)
        var_idx = min(var_idx, len(gene_var_list) - 1)
        while var_idx >= 0:
            _, var_id = gene_var_list[var_idx]
            if var_id.startswith("nv") \
                    or var_id in ht \
                    or var_id not in Links:
                var_idx -= 1
                continue
            if var_id in gene_var_maxrights \
                    and gene_var_maxrights[var_id] < left:
                break
            var_type, var_left, var_data = gene_vars[var_id]
            var_right = var_left
            if var_type == "deletion":
                var_right = var_left + int(var_data) - 1
            if (var_left >= left and var_left <= right) \
                    or (var_right >= left and var_right <= right):
                tmp_alleles |= set(Links[var_id])
            var_idx -= 1
        alleles -= tmp_alleles
        alleles &= set(count_per_read.keys())

        for allele in alleles:
            count_per_read[allele] += add

        return len(alleles)

    def add_stat(Gene_cmpt,
                 Gene_counts,
                 Gene_count_per_read,
                 include_alleles = set()):
        """
        Find the max of count in Gene_count_per_read
        to Gene_count

        And update Gene_cmpt by max allele-name then return the name
        """
        if len(Gene_count_per_read) <= 0:
            return ""
        max_count = max(Gene_count_per_read.values())
        cur_cmpt = set()
        for allele, count in Gene_count_per_read.items():
            if count < max_count:
                continue
            if len(include_alleles) > 0 \
                    and allele not in include_alleles:
                continue

            cur_cmpt.add(allele)
            if allele not in Gene_counts:
                Gene_counts[allele] = 1
            else:
                Gene_counts[allele] += 1

        if len(cur_cmpt) == 0:
            return ""

        verbose = 3
        if verbose >= 2:
            alleles = ["", ""]
            allele1_found = False
            allele2_found = False
            if alleles[0] != "":
                for allele, count in Gene_count_per_read.items():
                    if count < max_count:
                        continue
                    if allele == alleles[0]:
                        allele1_found = True
                    elif allele == alleles[1]:
                        allele2_found = True
                if allele1_found != allele2_found:
                    print((alleles[0],
                           Gene_count_per_read[alleles[0]]),
                          file=sys.stderr)
                    print((alleles[1],
                           Gene_count_per_read[alleles[1]]),
                          file=sys.stderr)
                    if allele1_found:
                        print("%s\tread_id %s - %d vs. %d]"
                               % (alleles[0],
                                  prev_read_id,
                                  max_count,
                                  Gene_count_per_read[alleles[1]]),
                              file=sys.stderr)
                    else:
                        print("%s\tread_id %s - %d vs. %d]"
                               % (alleles[1],
                                  prev_read_id,
                                  max_count,
                                  Gene_count_per_read[alleles[0]]),
                              file=sys.stderr)

        cur_cmpt = sorted(list(cur_cmpt))
        cur_cmpt = '-'.join(cur_cmpt)
        if not cur_cmpt in Gene_cmpt:
            Gene_cmpt[cur_cmpt] = 1
        else:
            Gene_cmpt[cur_cmpt] += 1

        # KIR2DP1*002-KIR2DP1*003-KIR2DP1*004-KIR2DP1*007-KIR2DP1*008-KIR2DP1*009-KIR2DP1*010
        return cur_cmpt

    def cmp_to_alt(cmp_list):
        # Remove mismatches due to unknown or novel variants
        cmp_list2 = []
        for cmp in cmp_list:
            cmp = deepcopy(cmp)
            type, pos, length = cmp[:3]
            if type == "match":
                if len(cmp_list2) > 0 and cmp_list2[-1][0] == "match":
                    cmp_list2[-1][2] += length
                else:
                    cmp_list2.append(cmp)
            elif type == "mismatch" and \
                    (cmp[3] == "unknown" or cmp[3].startswith("nv")):
                if len(cmp_list2) > 0 and cmp_list2[-1][0] == "match":
                    cmp_list2[-1][2] += 1
                else:
                    cmp_list2.append(["match", pos, 1])
            else:
                cmp_list2.append(cmp)

        # print(cmp_list)
        # print(cmp_list2)
        # [['match', 0, 171], ['mismatch', 171, 1, 'hv22'], ['match', 172, 115], ['mismatch', 287, 1, 'nv9'], ['match', 288, 6]]
        # [['match', 0, 171], ['mismatch', 171, 1, 'hv22'], ['match', 172, 122]]
        verbose = 2
        cmp_list_left, cmp_list_right, \
            cmp_left_alts, cmp_right_alts \
            = typing_common.identify_ambigious_diffs(ref_seq,
                                                     gene_vars,
                                                     Alts_left,
                                                     Alts_right,
                                                     Alts_left_list,
                                                     Alts_right_list,
                                                     cmp_list2,
                                                     verbose,
                                                     debug=False)

        mid_ht = []
        for cmp in cmp_list2[cmp_list_left:cmp_list_right + 1]:
            type = cmp[0]
            if type not in ["mismatch", "deletion", "insertion"]:
                continue
            var_id = cmp[3]
            mid_ht.append(var_id)

        # print(cmp_list, "->", cmp_list_left, cmp_list_right)
        # print("l", cmp_left_alts, "r", cmp_right_alts, "c", mid_ht)
        # [['match', 3512, 96], ['deletion', 3608, 1, 'hv1723'], ['match', 3609, 4], ['mismatch', 3613, 1, 'hv1725'], ['match', 3614, 1], ['mismatch', 3615, 1, 'hv1726'], ['mismatch', 3616, 1, 'hv1727'], ['mismatch', 3617, 1, 'hv1728'], ['mismatch', 3618, 1, 'hv1729'], ['mismatch', 3619, 1, 'hv1730'], ['match', 3620, 1], ['mismatch', 3621, 1, 'hv1732'], ['mismatch', 3622, 1, 'hv1733'], ['mismatch', 3623, 1, 'hv1736'], ['match', 3624, 1], ['mismatch', 3625, 1, 'hv1738'], ['match', 3626, 2], ['mismatch', 3628, 1, 'hv1739'], ['mismatch', 3629, 1, 'hv1740'], ['mismatch', 3630, 1, 'hv1742'], ['mismatch', 3631, 1, 'hv1743'], ['mismatch', 3632, 1, 'hv1744'], ['mismatch', 3633, 1, 'hv1745'], ['match', 3634, 1], ['mismatch', 3635, 1, 'hv1746'], ['match', 3636, 1], ['mismatch', 3637, 1, 'hv1747'], ['match', 3638, 15], ['mismatch', 3653, 1, 'hv1757'], ['match', 3654, 130], ['deletion', 3784, 1, 'hv1827'], ['match', 3785, 30]] -> 0 30
        # l ['3512'] r ['hv1841-hv1842-hv1843-hv1844-hv1845-hv1846-hv1847-hv1848-hv1850-hv1851-hv1852-hv1853-hv1854-hv1855-3816', '3814'] c ['hv1723', 'hv1725', 'hv1726', 'hv1727', 'hv1728', 'hv1729', 'hv1730', 'hv1732', 'hv1733', 'hv1736', 'hv1738', 'hv1739', 'hv1740', 'hv1742', 'hv1743', 'hv1744', 'hv1745', 'hv1746', 'hv1747', 'hv1757', 'hv1827']
        return cmp_left_alts, mid_ht, cmp_right_alts

    def cmp_to_hts(cmp_list):
        positive_hts   = set()
        cmp_left_alts, mid_ht, cmp_right_alts = cmp_to_alt(cmp_list)
        for left_id in range(len(cmp_left_alts)):
            left_ht = cmp_left_alts[left_id].split('-')
            left_ht += mid_ht
            for r in range(len(cmp_right_alts)):
                right_ht = cmp_right_alts[r].split('-')
                ht = left_ht + right_ht
                if len(ht) <= 0:
                    continue
                ht_str = '-'.join(ht)
                positive_hts.add(ht_str)
        return positive_hts

    # ----------------------
    # --- main ---
    # ----------------------
    for i in set(map(lambda i: i.split("*")[0], allele_vars.keys())):
        print(f"@RG\tID:{i}", file=bam_group_f)

    Gene_counts    = defaultdict(int)
    Gene_cmpt      = {}
    allele_stat    = defaultdict(int)
    all_reads      = []

    print("Start reading each record in Bam")

    for line_l, line_r, cmp_list_l, cmp_list_r in \
        readPair(alignment_fname, gene_vars, gene_var_list, ref_locus, ref_seq, base_locus):
        print(line_l, file=bam_group_f)
        print(line_r, file=bam_group_f)

        # ------------------------------
        # --- Processing Paired-read ---
        # ------------------------------
        left_positive_hts  = cmp_to_hts(cmp_list_l)
        right_positive_hts = cmp_to_hts(cmp_list_r)

        # linnil1 hisat-kir methods
        allele_count_p = defaultdict(int)
        allele_count_n = defaultdict(int)
        for positive_ht in left_positive_hts | right_positive_hts:
            count_allele_pn(positive_ht)
        allele_count_p_sort = sorted(allele_count_p.items(), key=lambda i: -i[1])
        allele_count_n_sort = sorted(allele_count_n.items(), key=lambda i: -i[1])
        all_reads.append({'p': allele_count_p_sort,
                          'n': allele_count_n_sort})
        # linnil1 method Done

        # hisat method
        Gene_count_per_read = {}
        # Gene_primary_exons_count_per_read = {}
        # Gene_exons_count_per_read         = {}
        for allele in gene_names:
            if allele.find("BACKBONE") != -1:
                continue
            if base_fname == "genome" and allele.find("GRCh38") != -1:
                continue
            # if allele in primary_exon_allele_rep_set:
            #     Gene_primary_exons_count_per_read[allele] = 0
            # if allele in allele_rep_set:
            #     Gene_exons_count_per_read[allele] = 0
            Gene_count_per_read[allele] = 0
        for positive_ht in left_positive_hts | right_positive_hts:
            # primary_exon_hts \
            #     = get_exon_haplotypes(positive_ht,
            #                           ref_primary_exons)
            # for exon_ht in primary_exon_hts:
            #     add_count(Gene_primary_exons_count_per_read,
            #               exon_ht,
            #               1)
            # exon_hts = get_exon_haplotypes(positive_ht,
            # for exon_ht in exon_hts:
            #     add_count(Gene_exons_count_per_read,
            #               exon_ht,
            #               1)
            num_allele = add_count(Gene_count_per_read,
                                   positive_ht,
                                   1)

        cur_cmpt = add_stat(Gene_cmpt,
                            Gene_counts,
                            Gene_count_per_read)
        if cur_cmpt:
            # print(cur_cmpt)
            allele_stat[len(cur_cmpt.split('-'))] += 1
            # if len(cur_cmpt.split('-')) == 1:
            #     print(line + "\tRG:Z:" + cur_cmpt, file=bam_group_f)
            #     print(prev_line + "\tRG:Z:" + cur_cmpt, file=bam_group_f)
        else:
            allele_stat[0] += 1
        # for read_id_, read_id_i, read_node in read_nodes:
        #     asm_graph.add_node(read_id_,
        #                        read_id_i,
        #                        read_node,
        #                        simulation)
        # clear
        # read_nodes    = []
        # read_var_list = []

    # if not num_reads:
    #     raise Exception("Cannot read bam file")

    # print(sorted(allele_stat.items(), key=lambda i: i[1], reverse=True)[:10])

    # Counting abundance (Same in hisat2)
    # Gene_counts = {"KIR2DL1*001": 12}
    Gene_counts = [[allele, count] for allele, count in Gene_counts.items()]
    Gene_counts = sorted(Gene_counts, key=lambda x: x[1], reverse=True)
    print(gene)
    for count_i in range(len(Gene_counts)):
        count = Gene_counts[count_i]
        print("%d %s (count: %d)" % (count_i + 1, count[0], count[1]))
        print("%d %s (count: %d)" % (count_i + 1, count[0], count[1]), file=report_f)
        if count_i >= 50:
            break
    Gene_prob = typing_common.single_abundance(Gene_cmpt)
    print("\n")

    for prob_i in range(len(Gene_prob)):
        prob = Gene_prob[prob_i]
        # if prob[1] < 0.01:
        #     break
        print("%d ranked %s (abundance: %.2f%%)" % (prob_i + 1, prob[0] , prob[1] * 100.0))
        print("%d ranked %s (abundance: %.2f%%)" % (prob_i + 1, prob[0] , prob[1] * 100.0), file=report_f)
        if prob_i >= 50:
            break

    # linnli1 method save tmp results
    json.dump(all_reads, open(typing_tmp_fname + f".{gene}.json", "w"))


def typingAllGene(alignment_bam_fname):
    genes = ["KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5", "KIR2DP1",
             "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5",
             "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DP1", "KIR3DS1"]

    # open files
    global alignment_fname, typing_tmp_fname, bam_group_f, report_f
    alignment_fname     = alignment_bam_fname
    typing_tmp_fname    = alignment_fname[:-4] + ".tmp"
    bam_group_fname     = alignment_fname[:-4] + ".tmp.sam"
    report_fname        = alignment_fname[:-4] + ".report"
    os.system(f"samtools view -H {alignment_fname} > {bam_group_fname}")
    report_f            = open(report_fname, "w")
    bam_group_f         = open(bam_group_fname, "a")

    # typing for each gene in hisat2
    print("Typing for", alignment_fname)
    for gene in refGenes:
        hisatTyping(gene)

    # close files
    bam_group_f.close()
    report_f.close()
    print("Report", report_fname)
    print("Postive Negative allele data", typing_tmp_fname)
    print("Postive Negative allele bam", bam_group_fname)


if __name__ == "__main__":
    # full_gg_path        = "./kir_split"
    # full_gg_path        = "./kir_merge"
    # alignment_fname     = "data/synSeq.hiseq.dp50.rl150.1.pair.split.bam"
    # python3 pipeline_typing.py ./kir_split data/synSeq.hiseq.dp50.rl150.1.merge.pair.bam
    full_gg_path        = sys.argv[1]
    alignment_fname     = sys.argv[2]
    if not os.path.exists(alignment_fname):
        print("Fail to read ", alignment_fname)
        sys.exit()
    hisatdataInit(full_gg_path)
    # typingAllGene(alignment_fname)
