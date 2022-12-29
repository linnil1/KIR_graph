from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
from graphkir.utils import readFromMSAs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def calcIntronExonPosition(seq: str, exon_regions: list[tuple[int, int]] = []) -> tuple[set[int], set[int]]:
    alls = set(range(len(seq)))
    gaps = set(i for i, base in enumerate(seq) if base == "-")
    exons = set()
    for left, right in exon_regions:
        for i in range(left, right):
            exons.add(i)
    exons -= gaps
    introns = alls - gaps - exons
    return exons, introns


def addNovelOnSeq(seq: str, want_novel_position: dict[str, int], exon_regions: list[tuple[int, int]] = []) -> tuple[str, list[tuple[int, str, str]]]:
    # position candidate
    exons_set, introns_set = calcIntronExonPosition(seq, exon_regions)
    alls = list(exons_set | introns_set)
    exons = list(exons_set)
    introns = list(introns_set)

    # select position
    positions: list[int] = []
    for want_region, want_variants_num in want_novel_position.items():
        if want_region == "all":
            positions.extend(np.random.choice(alls, want_variants_num, replace=False))
        elif want_region == "exon":
            positions.extend(np.random.choice(exons, want_variants_num, replace=False))
        elif want_region == "intron":
            positions.extend(np.random.choice(introns, want_variants_num, replace=False))
        else:
            raise NotImplementedError
    positions = sorted(positions)

    # apply
    apply_variants = []
    for pos in positions:
        base = np.random.choice(list(set("ACGT") - set(seq[pos])))
        apply_variants.append((pos, seq[pos], base))
        seq = seq[:pos] + base + seq[pos+1:]

    # return
    return seq, apply_variants


def addNovel(sequences: list[tuple[str, str]], exon_regions=dict[str, list[tuple[int, int]]], seed: int = 123) -> list[SeqRecord]:
    # read regions
    np.random.seed(seed)
    want_novel_position = {
        "exon": 2,
        "intron": 2,
    }
    new_sequences = []
    for name, seq in sequences:
        # TODO: double check if the new sequence
        seq, apply_variants = addNovelOnSeq(seq, want_novel_position, exon_regions[name])
        new_sequences.append(SeqRecord(
            Seq(seq.replace("-", "")),
            id=name + "".join(f"-{pos}{b}" for pos, a, b in apply_variants),
            description=",".join(f"{name}:{pos}{a}>{b}" for pos, a, b in apply_variants),
        ))
        print(new_sequences[-1].id, new_sequences[-1].description)
    return new_sequences


def addNovelFromFasta(input_fasta, input_bed, output_fasta, seed: int = 123):
    exon_regions = defaultdict(list)
    for line in open(input_bed):
        ref, pos1, pos2 = line.split("\t")
        exon_regions[ref].append((int(pos1), int(pos2)))
    sequences = list((record.id, str(record.seq)) for record in SeqIO.parse(input_fasta, "fasta"))
    new_sequences = addNovel(sequences, exon_regions, seed=seed)
    SeqIO.write(new_sequences, output_fasta, "fasta")


def addNovelFromMsa(input_fasta, msa_name, output_fasta, seed=123):
    genes = readFromMSAs(str(msa_name))
    name_to_gene = {}
    for gene, msa in genes.items():
        for allele_name in msa.alleles.keys():
            name_to_gene[allele_name] = gene

    answer_alleles = list(record.id.split("-")[0] for record in SeqIO.parse(input_fasta, "fasta"))
    sequences = []
    exon_regions = defaultdict(list)
    for allele_name in answer_alleles:
        msa = genes[name_to_gene[allele_name]]
        sequences.append((allele_name, msa.get(allele_name)))
        exons = [i.name for i in msa.blocks if i.type == "exon"]
        exon_regions[allele_name].extend(msa.get_block_interval(exon) for exon in exons)
    new_sequences = addNovel(sequences, exon_regions, seed=seed)
    SeqIO.write(new_sequences, output_fasta, "fasta")


def addNovelFromFastaWrap(input_name, bed_name, seed=123):
    output_name = input_name + f".noveli2e2_s{seed}"
    if Path(output_name + ".fa").exists():
        return output_name
    addNovelFromFasta(input_name + ".fa",
                      bed_name.format(input_name.template_args[0]) + ".bed",
                      output_name + ".fa",
                      seed=seed)
    return output_name


def addNovelFromMsaWrap(input_name, msa_name, seed=123):
    output_name = input_name + f".novelmsa_i2e2_s{seed}"
    if Path(output_name + ".fa").exists():
        return output_name
    addNovelFromMsa(input_name + ".fa",
                    msa_name,
                    output_name + ".fa",
                    seed=seed)
    return output_name


def updateNovelAnswer(input_name, old_name):
    output_path = input_name.replace(".{}", "_summary")
    answer_file = old_name.get_input_names([-1])[0].replace_wildcard("_summary")
    if Path(output_path + ".fa").exists():
        return input_nmae

    result_df = pd.read_csv(answer_file + ".csv", sep="\t", dtype="str")
    for name in input_name.get_input_names():
        id = name.template_args[-1]
        allele_names = list(i.id for i in SeqIO.parse(name + ".fa", "fasta"))
        result_df.loc[result_df['id'] == id, 'alleles'] = "_".join(allele_names)
    print(result_df)
    result_df.to_csv(output_path + ".csv", index=False, sep="\t")
    return input_name


if __name__ == "__main__":
    addNovelFromFasta("linnil1_syn/linnil1_syn_s44.00.fa", "linnil1_syn/linnil1_syn_s44.00.exon_region.bed", "test/add_novel.fasta")
