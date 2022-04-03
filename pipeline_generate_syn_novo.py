from collections import defaultdict
import os
import random
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyHLAMSA import Genemsa


def readMSA():
    index = "kir_merge_full"
    msa = Genemsa.load_msa(f"{index}.save.fa", f"{index}.save.gff")
    msa.gene_name = "KIR"
    return msa


def getPos(msa):
    gen_pos = msa._calculate_position()
    exon_pos   = list(zip([i[1] for i in msa.labels if i[0] == 'exon'], gen_pos[:-1], gen_pos[1:]))
    intron_pos = list(zip([i[1] for i in msa.labels if i[0] != 'exon'], gen_pos[:-1], gen_pos[1:]))
    return exon_pos, intron_pos


def getGene(msa):
    return sorted(set(map(lambda i: i.split("*")[0], msa.alleles.key())))


def getAllele(msa, gene):
    allele_names = filter(lambda i: not i.endswith("KIR*BACKBONE"), msa.alleles.keys())
    allele_names = filter(lambda i: i.startswith(gene + "*"), allele_names)
    return list(allele_names)


def setDenovo(msa, gene=None, allele=None, region="intron", pos=None):
    exon_pos, intron_pos = getPos(msa)

    if gene is None:
        gene = random.choice(getGene(msa))

    # allele name
    if allele is None:
        allele = random.choice(getAllele(msa, gene))

    while True:
        # intron pos
        if region == "intron":
            block = random.choice(intron_pos)
        elif region == "exon":
            block = random.choice(intron_pos)
        else:
            raise "region should be intron or exon"

        if pos is None:
            pos = random.randint(block[1], block[2])

        # replace
        if msa.alleles[allele][pos] == '-':
            pos = None
            continue
        new_seq = msa.alleles[allele]
        ori_nt = new_seq[pos]
        nt = random.choice(list(set("ATCG") - set(ori_nt)))
        new_seq = f"{new_seq[:pos]}{nt}{new_seq[pos]}"

        # check no same sequence
        if new_seq in msa.alleles.values():
            continue

        new_name = f"{allele}:{pos}:{ori_nt}>{nt}"
        print(new_name)
        return new_name, new_seq


def generateFastq(name):
    cmd = f"""
        PING/ping_paper_scripts/art_bin_MountRainier/art_illumina \
        -ss HS25 -na \
        -i {name}.fa \
        -l 150 \
        -f 50  \
        -m 400 \
        -s 10  \
        -sam -M \
        -rs 444 \
        -o {name}.read.
    """
    print(cmd)
    os.system(cmd)


def writeFasta(msa, allele_list, output_name):
    real_seqs = []
    allele_cn = defaultdict(int)
    for allele_name in allele_list:
        if not allele_cn[allele_name]:
            name_dup = allele_name
        else:
            name_dup = f"{allele_name}-{allele_cn[allele_name]}"
        real_seqs.append(
            SeqRecord(Seq(msa.alleles[allele_name].replace('-', '')),
                      id=name_dup,
                      description="")
        )
        allele_cn[allele_name] += 1
    SeqIO.write(real_seqs, f"{output_name}.fa", "fasta")


def writeAll(name, allele_list):
    summary.append({
        'id': name,
        'allele': '_'.join(sorted(allele_list)),
    })
    writeFasta(msa, allele_list, f"{folder}/{name}")
    generateFastq(f"{folder}/{name}")



# init
output_name = "linnil1_syn_novo"  # m = 400, N=100, seed=44
folder = output_name + "/" 
os.makedirs(folder, exist_ok=True)
random.seed(44)
summary = []
msa = readMSA()

# sample
name = output_name + ".1"
framework_1 = random.choice(getAllele(msa, "KIR3DL3"))
framework_2 = random.choice(getAllele(msa, "KIR3DL3"))
new_allele, new_seq = setDenovo(msa, gene="KIR2DL3")
msa.alleles[new_allele] = new_seq
normal_allele = random.choice(getAllele(msa, "KIR2DL3"))
alleles = [framework_1, framework_2, new_allele, normal_allele]
writeAll(name, alleles)

# sample2
name = output_name + ".2"
alleles = [framework_1, framework_2, new_allele, new_allele]
writeAll(name, alleles)

# sample3
name = output_name + ".3"
new_allele_exon, new_seq_exon = setDenovo(msa, region="exon", gene="KIR2DL3")
msa.alleles[new_allele_exon] = new_seq_exon
alleles = [framework_1, framework_2, new_allele_exon, normal_allele]
writeAll(name, alleles)

# summary
summary = pd.DataFrame(summary)
summary.to_csv(f"{folder}/summary.csv", sep="\t", header=False, index=False)
