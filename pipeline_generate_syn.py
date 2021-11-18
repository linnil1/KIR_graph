import pandas as pd
import random
import os
from Bio import SeqIO
from collections import defaultdict
import copy

# config
N = 10
output_name = "linnil1_syn"
output_name = "linnil1_syn_full"
folder = output_name + "/"
file_sequence = "kir_merge_sequences.fa"
file_sequence = "kir_merge_full_sequences.fa"

# read haplo
data = pd.read_csv("PING/ping_paper_scripts/synthetic_sequence/KIR_gene_haplotypes.csv")
data.index = data['hapID']

# read sequences
seqs = SeqIO.to_dict(SeqIO.parse(file_sequence, "fasta"))
gene_map_allele = defaultdict(list)
for name in seqs:
    # allele group
    if len(seqs[name]) > 4200:  # remove exon-only sequence, min 3DP1 4240
        # Note 2DL5A 2DLB are treat as 2DL5
        # gene_map_allele[name.split("*")[0]].append(name)
        gene_map_allele[name[:7]].append(name)
    else:
        print(f"remove {name} (len={len(seqs[name])})")

# run random
random.seed(444)
ans = []
for i in range(N):
    haplo = random.sample(list(data.index), k=2)
    alleles = []
    print(haplo)
    print(data.loc[haplo, :])

    for name, num in data.loc[haplo, :].iloc[:, 1:].sum().items():
        print(name, num)
        alleles.extend(random.choices(gene_map_allele[name], k=num))

    id = f"{output_name}.{i:02d}"
    ans.append({
        'id': id,
        'haplo': sorted(haplo),
        'allele': sorted(alleles),
    })

# save fasta
os.makedirs(folder, exist_ok=True)
for a in ans:
    real_seqs = []
    allele_cn = defaultdict(int)
    for allele_name in a['allele']:
        allele_cn[allele_name] += 1
        seq = copy.deepcopy(seqs[allele_name])
        seq.id = f"{seq.id}-{allele_cn[allele_name]}"
        seq.description = ""
        real_seqs.append(seq)

    SeqIO.write(real_seqs, f"{folder}/{a['id']}.fa", "fasta")

    # '-ef',  # error free
    """
    -l length
    -f depth
    -m mean size
    -s lenght devation
    -rs randomseed
    -sam -M  Output sam
    """
    cmd = f"""\
    PING/ping_paper_scripts/art_bin_MountRainier/art_illumina \
    -ss HS25 -na \
    -i {folder}/{a['id']}.fa \
    -l 150 \
    -f 50  \
    -m 200 \
    -s 10  \
    -sam -M \
    -rs 444 \
    -o {folder}/{a['id']}.read."""

    print(cmd)
    os.system(cmd)

# save summary
ans = pd.DataFrame(ans)
ans['haplo'] = ans['haplo'].str.join("_")
ans['allele'] = ans['allele'].str.join("_")
ans.to_csv(f"{folder}/summary.csv", sep="\t", header=False, index=False)
