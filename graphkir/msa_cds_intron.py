"""
The fillMissingIntrons function enable to fill the intron of exon-only sequences
by consensus of similar full-length sequences
"""
from typing import Iterable

from pyhlamsa import KIRmsa, Genemsa

from .utils import getAlleleField, limitAlleleField


def removeExonIncompleteSeq(msa: Genemsa) -> Genemsa:
    """Remove the sequence when E in exon (We only expect E in intron)"""
    remove_names = set()
    for msa_part in msa.split_block():
        if msa_part.blocks[0].type == "exon":
            for name, seq in msa_part.alleles.items():
                if "E" in seq:
                    remove_names.add(name)
                    print(f"remove {name} bcz E in {msa_part.blocks[0]}")
    print("remove", remove_names)
    msa = msa.remove_allele(remove_names)
    return msa


def searchNearestName(full_names: Iterable[str], target_name: str) -> list[str]:
    """
    Find the alleles that simiar to the target_allele

    If return empty list, it indicate they doesn't exist any allele
    that has at least three digits similar to target
    """
    if not target_name[-1].isdigit():
        target_name = target_name[:-1]
    field = len(getAlleleField(target_name))
    assert field in [3, 5, 7]

    while True:
        nearest_name = [name for name in full_names if name.startswith(target_name)]
        if nearest_name:
            return nearest_name
        # else -> cannot find
        if field == 3:
            break
        # 5 or 7
        field -= 2
        target_name = limitAlleleField(target_name, field)

    return []


def getNearestConsensus(msa: Genemsa, target_names: list[str]) -> str:
    """Generate the consensus by the alleles in target_names"""
    if not target_names:
        # all sequence
        return msa.select_complete().get_consensus(include_gap=True)
    new_msa = msa.select_allele(target_names)
    return new_msa.get_consensus(include_gap=True)


def fillByConsensus(seq: str, consensus: str) -> str:
    """Fill the seq's E by consensus"""
    return "".join([seq[i] if seq[i] != "E" else consensus[i] for i in range(len(seq))])


def fillByNearestName(msa: Genemsa) -> Genemsa:
    """Fill the exon-only sequence by the most similar alleles"""
    new_msa = msa.copy(copy_allele=False)
    full_names = msa.select_complete().alleles.keys()
    exon_names = msa.select_incomplete().alleles.keys()

    for name in full_names:
        new_msa.append(name, msa.get(name))

    for name in exon_names:
        nearest_names = searchNearestName(full_names, name)
        print(name, nearest_names)
        consensus = getNearestConsensus(msa, nearest_names)
        seq = fillByConsensus(msa.get(name), consensus)
        new_msa.append(name + "e", seq)

    return new_msa


def fillMissingIntrons(genes: dict[str, Genemsa]) -> dict[str, Genemsa]:
    """Fill the exon-only gap by its similar alleles for all genes"""
    new_kir = {}
    for gene, msa in genes.items():
        print(gene)
        msa = removeExonIncompleteSeq(msa)
        msa = fillByNearestName(msa)
        new_kir[gene] = msa
    return new_kir


if __name__ == "__main__":
    from graphkir.kir_msa import saveAllMsa

    kir = KIRmsa(filetype=["nuc", "gen"], version="2100")
    genes = fillMissingIntrons(kir.genes)
    saveAllMsa(genes, "index5/kir_2100_withexon")
