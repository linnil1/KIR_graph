import pandas as pd
import random
import os
from Bio import SeqIO
from collections import defaultdict

# read haplo
data = pd.read_csv("PING/ping_paper_scripts/synthetic_sequence/KIR_gene_haplotypes.csv")
data.index = data['hapID']

# read sequences
seqs = SeqIO.to_dict(SeqIO.parse("kir_merge_sequences.fa", "fasta"))
gene_map_allele = defaultdict(list)
for name in seqs:
    # gene_map_allele[name.split("*")[0]].append(name)
    # allele group
    gene_map_allele[name[:7]].append(name)

# config
N = 10
folder = "linnil1_syn1/"

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

    id = f"linnil1_syn.{i:02d}"
    ans.append({
        'id': id,
        'haplo': sorted(haplo),
        'allele': sorted(alleles),
    })

# save fasta
os.makedirs(folder, exist_ok=True)
for a in ans:
    SeqIO.write([seqs[i] for i in a['allele']], f"{folder}/{a['id']}.fa", "fasta")

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
