"""
* read DB -> MSA
* split/merge the MSA
* realign the MSA when merge (muscle, clustalo)
"""
from glob import glob
from typing import Callable
from itertools import chain
from collections import defaultdict

from Bio import SeqIO, SeqRecord, Align
from pyhlamsa import msaio, Genemsa, KIRmsa
from .utils import runDocker
from .msa_cds_intron import fillMissingIntrons


# TODO:
# * intron56
# * compare muscle and clustalo
GenesMsa = dict[str, Genemsa]
BlockMsa = dict[str, list[SeqRecord]]
BlockFile = dict[str, str]
kir_block_name = [
    "5UTR", "exon1", "intron1", "exon2", "intron2", "exon3", "intron3",
    "exon4", "intron4", "exon5", "intron5", "exon6", "intron6", "exon7",
    "intron7", "exon8", "intron8", "exon9", "3UTR"]


def saveAllMsa(genes: GenesMsa, prefix: str) -> None:
    """
    Save each gene's MSA

    Args:
        genes: A dictionary of gene's name and its MSA
        prefix: Save the MSA into filename with the prefix
    """
    for gene_name, msa in genes.items():
        msa = msa.shrink()
        msa.append(f"{gene_name}*BACKBONE",
                   msa.get_consensus(include_gap=False))
        msa.set_reference(f"{gene_name}*BACKBONE")
        msaio.to_bam(msa,   f"{prefix}.{gene_name}.bam")  # noqa: E201
        msaio.to_gff(msa,   f"{prefix}.{gene_name}.gff")
        msaio.save_msa(msa, f"{prefix}.{gene_name}.fa",
                            f"{prefix}.{gene_name}.json")


def readDB(full_length_only: bool = True, version: str = "latest") -> GenesMsa:
    """
    Read IPD-KIR database

    Args:
        version: The version of IPD-KIR (default: latest)
    Returns:
        genes: A dictionary of gene's name and its MSA
    """
    if full_length_only:
        kir = KIRmsa(filetype=["gen"], version=version)
    else:
        kir = KIRmsa(filetype=["gen", "nuc"], version=version)
    return kir.genes


def readFromMSAs(prefix: str) -> GenesMsa:
    """
    Read MSAs via `{prefix}.*.json`

    Returns:
        genes: A dictionary of gene's name and its MSA
    """
    genes = {}
    for filename in glob(prefix + ".*.json"):
        split_name = filename[len(prefix) + 1:].split('.')
        # prefix.anotherprefix.*.json will not included
        if len(split_name) != 2:
            continue
        print("read", filename)
        gene = split_name[0]
        genes[gene] = msaio.load_msa(filename[:-5] + ".fa", filename)
    return genes


def removeBackbone(genes: GenesMsa) -> GenesMsa:
    """ Remove {gene}*BACKBONE alleles """
    for gene, msa in genes.items():
        if f"{gene}*BACKBONE" in msa.alleles:
            msa.remove(f"{gene}*BACKBONE")
    return genes


def splitMsaToBlocks(genes: GenesMsa) -> BlockMsa:
    """
    Split MSA by blocks (i.e 5UTR, intron1, exon1, ...)

    Args:
        genes: A dictionary of gene's name and its MSA
    Returns:
        blocks: A dictionary of block's name and its sequences
    """
    blocks = defaultdict(list)
    for msa in genes.values():
        # TODO: intron3/4 -> intron4
        for msa_part in msa.split():
            block_name = msa_part.blocks[0].name
            if block_name == "intron3/4":
                block_name = "intron3"
            elif block_name == "intron5/6":
                block_name = "intron5"
            blocks[block_name].extend(msa_part.to_records(gap=False))
    assert set(kir_block_name) == blocks.keys()
    return blocks


def blockToFile(blocks: BlockMsa, tmp_prefix: str = "tmp") -> BlockFile:
    """
    Save blocks' sequences into fasta

    Args:
        blocks: A dictionary of block's name and its sequences
    Returns:
        files: A dictionary of block's name and its filename
    """
    files = {}
    for block_name, recs in blocks.items():
        filename = tmp_prefix + "." + block_name
        files[block_name] = filename
        recs = [rec for rec in recs if len(rec.seq)]  # remove 0 length sequences
        SeqIO.write(recs, filename + ".fa", "fasta")
    return files


def realignBlock(files: BlockFile, method: str = "clustalo", threads: int = 1) -> BlockFile:
    """
    Run MSA tools for the files in `files`

    Args:
        files: A dictionary of block's name and its filename
    Returns:
        files: A dictionary of block's name and its filename after realign
    """
    files_aligned = {}
    for block_name, f in files.items():
        if method == "clustalo":
            result = clustalo(f, threads)
        elif method == "muscle":
            result = muscle(f, threads)
        else:
            raise NotImplementedError
        files_aligned[block_name] = result
    return files_aligned


def fileToBlock(files: BlockFile) -> BlockMsa:
    """
    Load blocks' fasta into sequences

    Args:
        files: A dictionary of block's name and its filename
    Returns:
        blocks: A dictionary of block's name and its sequences
    """
    return {
        block_name: list(SeqIO.parse(f + ".fa", "fasta"))
        for block_name, f in files.items()
    }


def mergeBlockToMsa(blocks: BlockMsa) -> Genemsa:
    """
    Merge blocks' sequences into MSA

    The sequences should has

    * same order with kir_block_name
    * which Genemsa.assume_label will works

    Args:
        blocks: A dictionary of block's name and its sequences
    Returns:
        Genemsa: An aligned MSA
    """
    # concat together
    msa = None
    current_alleles: set[str] = set()
    for block_name in kir_block_name:
        # seqrecord -> Genemsa
        blk = Genemsa.from_MultipleSeqAlignment(
                Align.MultipleSeqAlignment(blocks[block_name]))
        # first
        if msa is None:
            msa = blk
            current_alleles = set(msa.alleles.keys())
            continue

        # fill miss alleles
        for allele in blk.alleles.keys() - current_alleles:
            msa.append(allele, '-' * msa.get_length())
        for allele in current_alleles - blk.alleles.keys():
            blk.append(allele, '-' * blk.get_length())

        # merge
        msa += blk
        current_alleles = msa.alleles.keys()

    assert msa
    msa = msa.assume_label("gen")
    return msa


def isEqualMsa(msas: GenesMsa, msa: Genemsa) -> bool:
    assert set(msa.alleles.keys()) \
        == set(chain.from_iterable(msa_old.alleles.keys() for msa_old in msas.values()))
    for msa_old in msas.values():
        for name, seq in msa_old.alleles.items():
            if seq.replace("-", "") != msa.get(name).replace("-", ""):
                print(name)
                for i, j in zip(seq.replace("-", ""), msa.get(name).replace("-", "")):
                    print(i, j)
            assert seq.replace("-", "") == msa.get(name).replace("-", "")
    return True


def mergeMSA(genes: GenesMsa,
             method: str = "clustalo",
             tmp_prefix: str = "tmp",
             threads: int = 1) -> Genemsa:
    """
    Merge multiple MSA into one single MSA

    Args:
        genes: A dictionary of gene's name and its MSA
        method: The tools for realign the MSA Options: [clustalo, muscle]

    Returns:
        Genemsa: An aligned MSA
    """
    blocks = splitMsaToBlocks(genes)
    files = blockToFile(blocks, tmp_prefix=tmp_prefix)
    files = realignBlock(files, method, threads=threads)
    blocks = fileToBlock(files)
    msa = mergeBlockToMsa(blocks)
    isEqualMsa(genes, msa)
    return msa


def muscle(name: str, threads: int = 1) -> str:
    """
    Construct MSA with muscle
    (Input and output are filename without suffix)
    """
    runDocker("muscle",
              f"muscle -align {name}.fa -threads {threads}"
              f"       -output {name}.muscle.fa")
    return name + ".muscle"


def clustalo(name: str, threads: int = 1) -> str:
    """
    Construct MSA with clustalo
    (Input and output are filename without suffix)
    """
    runDocker("clustalo",
              f"clustalo --infile {name}.fa -o {name}.clustalo.fa"
              f"         --outfmt fasta --threads {threads} --force")
    return name + ".clustalo"


def buildKirMsa(mode: str, prefix: str, version: str = "2100",
                input_msa_prefix: str = "",
                full_length_only: bool = True,
                mergeMSA: Callable = mergeMSA,
                threads: int = 1) -> None:
    """
    Read KIR from database and save MSA into files with prefix

    e.g. `{prefix}.KIR2DL1.json`, `{prefix}.KIR2DL1.fa`

    Args:
        mode:
            * split: 17 MSA
            * ab: 16 MSA (merge KIR2DL5A/B)
            * ab_2dl1s1: 15 genes (merge KIR2DL5A/B, KIR2DL1/S1)
            * merge: 1 MSA (merge all)

        prefix: The prefix of filename for saving
        version: The version of IPD-KIR database
        input_msa_prefix: The prefix of MSA, we will read MSA from the path
                          instead of downloing `version` DB
    """
    if input_msa_prefix:
        genes = readFromMSAs(input_msa_prefix)
        genes = removeBackbone(genes)
    else:
        genes = readDB(full_length_only=full_length_only, version=version)
        if not full_length_only:
            genes = fillMissingIntrons(genes)

    if mode == "split":
        genes['KIR2DL5A'] = genes['KIR2DL5'].select_allele("KIR2DL5A.*")
        genes['KIR2DL5B'] = genes['KIR2DL5'].select_allele("KIR2DL5B.*")
        del genes['KIR2DL5']
    elif mode == "ab":
        pass
    elif mode == "merge":
        genes = {"KIR": mergeMSA(genes,
                                 method="clustalo",
                                 tmp_prefix=prefix + ".tmp",
                                 threads=threads)}
    elif mode == "ab_2dl1s1":
        genes_for_merge = {
            "KIR2DL1": genes.pop("KIR2DL1"),
            "KIR2DS1": genes.pop("KIR2DS1"),
        }
        genes["KIR2DL1S1"] = mergeMSA(genes_for_merge,
                                      method="muscle",
                                      tmp_prefix=prefix + ".tmp")
    else:
        raise NotImplementedError

    saveAllMsa(genes, prefix)
