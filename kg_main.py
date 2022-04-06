import uuid
import os
from glob import glob
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

from pyHLAMSA import KIRmsa, Genemsa
import pandas as pd
from Bio import SeqIO, AlignIO
import kg_build_index


thread = 30
kir_block_name = ["5UTR", "exon1", "intron1", "exon2", "intron2", "exon3",
                  "intron3", "exon4", "intron4", "exon5", "intron5", "exon6",
                  "intron6", "exon7", "intron7", "exon8", "intron8", "exon9",
                  "3UTR"]


def run_dk(image, cmd):
    """ run docker container """
    name = str(uuid.uuid4()).split("-")[0]
    run("docker run -it --rm --security-opt label=disable -u root -w /app "
        f"--name {name} -v $PWD:/app {image} {cmd}")


def run(cmd):
    """ wrap os.system """
    print(cmd)
    os.system(cmd)


def kir_to_multi_msa(split_2DL5=False):
    global index
    kir = KIRmsa(filetype=["gen"], version="2100")

    if split_2DL5:
        kir.genes['KIR2DL5A'] = kir.genes['KIR2DL5'].select_allele("KIR2DL5A.*")
        kir.genes['KIR2DL5B'] = kir.genes['KIR2DL5'].select_allele("KIR2DL5B.*")

    for gene_name, kir_msa_align in kir.genes.items():
        kir_msa_align = kir_msa_align.shrink()
        kir_msa_align.append(f"{gene_name}*BACKBONE",
                             kir_msa_align.get_consensus(include_gap=False))
        kir_msa_align.save_bam(f"{index}.{gene_name}.save.bam", ref_allele=f"{gene_name}*BACKBONE")
        kir_msa_align.save_gff(f"{index}.{gene_name}.save.gff", ref_allele=f"{gene_name}*BACKBONE")
        kir_msa_align.save_msa(f"{index}.{gene_name}.save.fa", f"{index}.{gene_name}.save.json")


def kir_to_single_msa(method="clustalo"):
    global index, suffix
    # Step1: break
    # kir_to_single_msa_break_block()
    suffix += ".tmp"

    # Step2: Build msa
    # kir_to_single_msa_realign(method=method)

    # Step3:
    msa = kir_to_single_msa_merge_block(ends_with=f"{method}.fa")

    # save the msa
    msa = msa.shrink().reset_index()
    msa.append(f"KIR*BACKBONE", msa.get_consensus(include_gap=False))
    msa.save_bam(f"{index}.save.bam", ref_allele=f"KIR*BACKBONE")
    msa.save_gff(f"{index}.save.gff", ref_allele=f"KIR*BACKBONE")
    msa.save_msa(f"{index}.save.fa", f"{index}.save.json")


def kir_to_single_msa_break_block():
    global index, suffix

    kir = KIRmsa(filetype=["gen"], version="2100")
    blocks = {block_name: [] for block_name in kir_block_name}
    block_length = []

    for gene_name, msa in kir.genes.items():
        b_length = {'name': gene_name}
        for i in range(len(msa.blocks)):
            if msa.blocks[i].name == "intron3/4":
                block_name = "intron3"
            elif msa.blocks[i].name == "intron5/6":
                block_name = "intron5"
            else:
                block_name = msa.blocks[i].name
            b_length[block_name] = msa.blocks[i].length
            blocks[block_name].extend(msa.select_block([i]).to_fasta(gap=False))
        block_length.append(b_length)

    print(pd.DataFrame(block_length))
    for block_name, seqs in blocks.items():
        SeqIO.write(filter(lambda i: len(i.seq), seqs),
                    open(f"{index}.tmp.{block_name}.fa", "w"),
                    "fasta")
    assert set(kir_block_name) == set(blocks.keys())


def kir_to_single_msa_realign(method):
    global index, suffix
 
    blocks = glob(f"{index}{suffix}.*.fa")
    with ProcessPoolExecutor(max_workers=thread) as executor:
        for name in blocks:
            name = os.path.splitext(name)[0]
            if method == "clustalo":
                executor.submit(clustalo, name)
            elif method == "muscle":
                executor.submit(muscle, name)
            else:
                raise NotImplementedError


def clustalo(name):
    run_dk("quay.io/biocontainers/clustalo:1.2.4--h1b792b2_4",
           f"clustalo --infile {name}.fa -o {name}.clustalo.fa "
           f"--outfmt fasta --threads 2 --force")

def muscle(name):
    run_dk("quay.io/biocontainers/muscle:5.1--h9f5acd7_1",
           f"muscle -align {name}.fa -threads 2"
           f"      -output {name}.muscle.fa")


def kir_to_single_msa_merge_block(ends_with):
    global index, suffix

    # Read all blocks
    blocks = glob(f"{index}{suffix}.*{ends_with}")
    block_msa = {}
    allele_names = set()
    for name in blocks:
        block_name = name[len(index) + len(suffix) + 1: -len(ends_with) - 1]
        block_msa[block_name] = Genemsa.from_MultipleSeqAlignment(
                AlignIO.read(name, "fasta"))
        allele_names.update(set(block_msa[block_name].get_sequence_names()))
    assert set(kir_block_name) == set(block_msa.keys())

    # add gap if the allele sequence in the block is empty
    for msa in block_msa.values():
        for an in allele_names:
            if an not in msa.alleles.keys():
                msa.append(an, '-' * msa.get_length())

    # concat together
    newmsa = Genemsa("KIR", "gen")
    newmsa.alleles = {name: "" for name in allele_names}
    for block_name in kir_block_name:
        msa = block_msa[block_name]
        """
        # not need to rename because assume_label is useful
        msa.blocks[0].name = block_name
        if "exon" in block_name:
            msa.blocks[0].type = "exon"
        elif block_name == "3UTR":
            msa.blocks[0].type = "three_prime_UTR"
        elif block_name == "5UTR":
            msa.blocks[0].type = "five_prime_UTR"
        """
        newmsa += msa
    newmsa = newmsa._assume_label().reset_index()

    # double check
    kir = KIRmsa(filetype=["gen"], version="2100")
    for msa in kir.genes.values():
        for name, seq in msa.alleles.items():
            assert seq.replace("-", "") == newmsa.get(name).replace("-", "")
    return newmsa


def samtobam():
    index_name = index.split("/")[-1]
    for name in samples:
        name += f".{index_name}{suffix}"
        run(f"samtools sort {name}.sam -o {name}.bam")
        run(f"samtools index {name}.bam")


def bamFilter(flag, new_suffix):
    index_name = index.split("/")[-1]
    for name in samples:
        name += f".{index_name}{suffix}"
        run(f"samtools view {name}.bam {flag} -h | samtools sort - -o {name}{new_suffix}.bam")
        run(f"samtools index {name}{new_suffix}.bam")


def hisatMap():
    index_name = index.split("/")[-1]
    for name in samples:
        name += f"{suffix}"
        f1, f2 = name + ".read1.fq", name + ".read2.fq"
        run(f""" \
            hisat2 --threads {thread} -x {index}.graph -1 {f1} -2 {f2} \
            --no-spliced-alignment --max-altstried 64 --haplotype \
            > {name}.{index_name}.sam
        """)
    samtobam()


samples = [f"data/linnil1_syn_wide.{i:02d}" for i in range(100)]
samples = samples[:1]

suffix = ""
# mkdir index
index = "index/kir_2100_raw"
# kir_to_multi_msa()
# kg_build_index.main(index)
index += ".mut01"
# hisatMap()
bamFilter("-f 0x2", ".nosingle")
# suffix += ".nosec"


# bamFilter("-f 0x2 -F 256", ".nosingle")



index = "index/kir_2100_ab"
# kir_to_multi_msa(split_2DL5=True)

index = "index/kir_2100_merge"
# kir_to_single_msa()

index = "index/kir_2100_muscle"
# kir_to_single_msa(method='muscle')
