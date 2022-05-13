import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO, AlignIO
from pyHLAMSA import KIRmsa, Genemsa

from kg_utils import runDocker, getSamples, threads


kir_block_name = [
    "5UTR", "exon1", "intron1", "exon2", "intron2", "exon3", "intron3",
    "exon4", "intron4", "exon5", "intron5", "exon6", "intron6", "exon7",
    "intron7", "exon8", "intron8", "exon9", "3UTR"]


def saveAllMsa(name, genes):
    """
    Args:
      index(str): Index name
      genes(dict[str, Genemsa]): A dictionary of gene_name and gene_msa
    """
    for gene_name, msa in genes.items():
        msa = msa.shrink()
        msa.append(f"{gene_name}*BACKBONE",
                   msa.get_consensus(include_gap=False))
        msa.set_reference(f"{gene_name}*BACKBONE")
        msa.save_bam(f"{name}.save.{gene_name}.bam")
        msa.save_gff(f"{name}.save.{gene_name}.gff")
        msa.save_msa(f"{name}.save.{gene_name}.fa",
                     f"{name}.save.{gene_name}.json")
    return ".save"


def kirToMultiMsa(index="index", split_2DL5=False):
    """
    Read KIR and output {index}.save.{gene_name}.xx
    """
    index = f"{index}/kir_2100_raw"
    if getSamples(index + ".save", strict=False):
        return index + ".save"

    kir = KIRmsa(filetype=["gen"], version="2100")
    if split_2DL5:
        kir.genes['KIR2DL5A'] = kir.genes['KIR2DL5'].select_allele("KIR2DL5A.*")
        kir.genes['KIR2DL5B'] = kir.genes['KIR2DL5'].select_allele("KIR2DL5B.*")
        del kir.genes['KIR2DL5']
        index = f"{index}/kir_2100_ab"

    index += saveAllMsa(index, kir.genes)
    return index


def kirToSingleMsa(index="index", method="clustalo"):
    if method == "clustalo":
        index = f"{index}/kir_2100_merge"
    elif method == "muscle":
        index = f"{index}/kir_2100_merge_muscle"
    else:
        raise NotImplementedError

    assign_intron_method = 0
    assign_intron_method = 1
    if assign_intron_method == 1:
        index += "_assign1"
    if getSamples(index + ".save", strict=False):
        return index + ".save"

    # Step1: break blocks
    kir = KIRmsa(filetype=["gen"], version="2100")
    if assign_intron_method == 0:
        blocks = genesToBlocks(kir.genes)
    elif assign_intron_method == 1:
        blocks = genesToBlocks(kir.genes, intron34="intron4", intron56="intron6")
    else:
        raise NotImplementedError

    # TODO: block-seq plot
    index += ".tmp"
    suffix = ""
    assert set(kir_block_name) == set(blocks.keys())
    for block_name, seqs in blocks.items():
        filename = f"{index}.{block_name}.fa"
        SeqIO.write(seqs, open(filename, "w"), "fasta")

    # Step2: Realign, build msa
    suffix += sequencesRealign(index, method=method)

    # Step3: double check and save the msa
    blocks = readBlocks(index, suffix=suffix)
    msa = blocksToMsa(blocks)
    kir = KIRmsa(filetype=["gen"], version="2100")  # test
    for msa_old in kir.genes.values():
        for name, seq in msa_old.alleles.items():
            assert seq.replace("-", "") == msa.get(name).replace("-", "")

    index = ".".join(index.split(".")[:-1])  # remove .tmp
    index += saveAllMsa(index[:-4], {'KIR': msa})
    return index


def genesToBlocks(genes, intron34="intron3", intron56="intron5"):
    """
    Args:
      genes(dict[str, Genemsa]): A dictionary of gene_name and gene_msa
    Returns:
      blocks(dict[str, list[str]]): A dictionary of block_name and list of sequences
    """
    blocks = defaultdict(list)
    for gene_name, msa in genes.items():
        for i in range(len(msa.blocks)):
            if msa.blocks[i].name == "intron3/4":
                block_name = intron34
            elif msa.blocks[i].name == "intron5/6":
                block_name = intron56
            else:
                block_name = msa.blocks[i].name
            blocks[block_name].extend(
                filter(lambda i: len(i.seq),
                       msa.select_block([i]).to_fasta(gap=False))
            )
    return blocks


def sequencesRealign(name, method):
    """
    Input:
      name(str): We will realign all {name}.*.fa with method
    """
    files = getSamples(name, ".fa")
    with ProcessPoolExecutor(max_workers=threads) as executor:
        for name in files:
            name = os.path.splitext(name)[0]
            if method == "clustalo":
                executor.submit(clustalo, name)
            elif method == "muscle":
                executor.submit(muscle, name)
            else:
                raise NotImplementedError
    return "." + method


def clustalo(name):
    runDocker("clustalo",
              f"clustalo --infile {name}.fa -o {name}.clustalo.fa "
              f"         --outfmt fasta --threads 4 --force")
    return ".clustalo"


def muscle(name):
    runDocker("muscle",
              f"muscle -align {name}.fa -threads 4 "
              f"       -output {name}.muscle.fa")
    return ".muscle"


def readBlocks(name, suffix):
    """
    Reads all {name}.{block}.{suffix}.fa

    Return:
      blocks(dict[str, Genemsa]): Dictionary of block_name and Genemsa
    """
    files = getSamples(name, suffix + ".fa", return_name=True)

    # Read all blocks
    block_msa = {}
    for filename, block_name in files:
        block_msa[block_name] = Genemsa.from_MultipleSeqAlignment(
            AlignIO.read(filename, "fasta"))

    # get all sequence names
    allele_names = set()
    for msa in block_msa.values():
        allele_names.update(msa.get_sequence_names())

    # add gap if the allele sequence in the block is empty
    for msa in block_msa.values():
        for an in allele_names:
            if an not in msa.alleles.keys():
                msa.append(an, '-' * msa.get_length())
    return block_msa


def blocksToMsa(blocks):
    """
    We assume the data
    * same order with kir_block_name
    * assume_label give same result of original kir block
    """
    # concat together
    newmsa = None
    for block_name in kir_block_name:
        if newmsa is None:
            newmsa = blocks[block_name]
        else:
            newmsa += blocks[block_name]
    newmsa = newmsa.assume_label("gen")
    return newmsa


def kirMerge2dl1s1(index="index"):
    """ Build msa from IPDKIR but merge 2dl1 and 2ds1 """
    kir2dls1 = ["KIR2DS1", "KIR2DL1"]
    index = f"{index}/kir_2100_2dl1s1"
    if getSamples(index + ".save", strict=False):
        return index + ".save"

    # Step1: Extract 2DL1 2DS1
    index += ".tmp"
    suffix = ""
    kir = KIRmsa(filetype=["gen"], version="2100")
    blocks = genesToBlocks({i: kir.genes[i] for i in kir2dls1})
    for block_name, seqs in blocks.items():
        filename = f"{index}.{block_name}.fa"
        SeqIO.write(seqs, open(filename, "w"), "fasta")

    # Step2: Build msa
    suffix += sequencesRealign(index, "muscle")
    blocks = readBlocks(index, suffix=suffix)
    msa = blocksToMsa(blocks)

    # Step3: double check and merge it
    for gene_name in kir2dls1:
        msa_old = kir.genes[gene_name]
        for name, seq in msa_old.alleles.items():
            assert seq.replace("-", "") == msa.get(name).replace("-", "")
        del kir.genes[gene_name]
    kir.genes["KIR2DL1S1"] = msa

    # Step4: save to msa
    index = index[:-4]
    index += saveAllMsa(index, kir.genes)
    return index


if __name__ == "__main__":
    index = "index1"
    # index += kirToMultiMsa(index, split_2DL5=False)
    index = kirToSingleMsa(index, method="muscle")
    # index = kirToSingleMsa(index, method="clustalo")
    # index = kirMerge2dl1s1(index)
    print(index)
