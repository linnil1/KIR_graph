import re
import os
from pyHLAMSA import Genemsa, BlockInfo
import kg_build_index
from kg_utils import runDocker
from kg_main import hisatMap
from kg_typping import readPair


def createGene():
    index = "index1/test_gen.save"
    if os.path.exists(index + ".fa"):
        return index
    test_seq0 = "TTTCCCTGGGAATTTAAATCATTTTAGCTGGTTCTGCTGTAATACTAGAAATACAAGCATGAAAAATTCTAATGGTTTATTAGTCACAATGACTCCGAAAACATTAATAATACCTATTAGATACTTTGCATATTACACAGGAAGAAGAGTTTGAATCTCAGATAAAAACAATAAAAATACATGAAAAGTCTTTCATGTTAGCACAGATTTTAGGCATCTCGTGTTCGGATAAAAATACATGAAAAGTCTTTCACGTTAGCACAGATTTTAGGCATCTTGTGTTCGGGAGGTTGGATCTGA"
    test_seq1 = "TTTCCCTGAGAATTTAAATCATTTTAGCTGGTTCTGCTGTAATACTAGAAATACAAGCATGAAAAATTCTAATGGTTTATTAGTCACAATGACTCCGAAAACATTAATAATACCTATTAGATACTTTGCATATTACACAGGAAGAAGAGTTTGAATCTCAGATAAAAACAATAAAAATACATGAAAAGTCTTTCATGTTAGCACAGATTTTAGGCATCTCGTGTTCGGATAAAAATACATGAAAAGTCTTTCACGTTAGCACAGATTTTAGGCATCTTGTGTTCGGGAGGTTGGATCTGA"
    test_seq2 = "TTTCCCTGGGAATATAAATCATTTTAGCTGGTTCTGCTGTAATACTAGAAATACAAGCATGAAAAATTCTAATGGTTTATTAGTCACAATGACTCCGAAAACATTAATAATACCTATTAGATACTTTGCATATTACACAGGAAGAAGAGTTTGAATCTCAGATAAAAACAATAAAAATACATGAAAAGTCTTTCATGTTAGCACAGATTTTAGGCATCTCGTGTTCGGATAAAAATACATGAAAAGTCTTTCACGTTAGCACAGATTTTAGGCATCTTGTGTTCGGGAGGTTGGATCTGA"
    test_seq3 = "TTTCCCTGAGAATATAAATCATTTTAGCTGGTTCTGCTGTAATACTAGAAATACAAGCATGAAAAATTCTAATGGTTTATTAGTCACAATGACTCCGAAAACATTAATAATACCTATTAGATACTTTGCATATTACACAGGAAGAAGAGTTTGAATCTCAGATAAAAACAATAAAAATACATGAAAAGTCTTTCATGTTAGCACAGATTTTAGGCATCTCGTGTTCGGATAAAAATACATGAAAAGTCTTTCACGTTAGCACAGATTTTAGGCATCTTGTGTTCGGGAGGTTGGATCTGA"
    msa = Genemsa("test")
    msa.alleles = {
        'test*BACKBONE': test_seq0,
        'test*001': test_seq1,
        'test*002': test_seq2,
    }
    msa.blocks = [BlockInfo(length=len(test_seq0), name="exon1", type="exon")]
    msa = msa.reset_index()
    print(msa.format_alignment_diff())
    msa.save_msa(index + ".test.fa", index + ".test.json")

    msa.append('test*003', test_seq3)
    msa.save_fasta(index + ".fa")
    return index


def createRead(index):
    if os.path.exists(index + ".read.1.fq"):
        return index
    runDocker("alpine", f"""\
        ./art_illumina \
        -ef \
        -ss HS25 \
        -i {index}.fa \
        -l 150 \
        -m 300 \
        -f 30 \
        -s 10 \
        -sam -na\
        -rs 444 \
        -o {index}.read. \
    """)
    return index


def inspectMapping(sample_index):
    for r1, r2 in readPair(sample_index + ".bam"):
        print(re.findall("Zs:Z:(.*) ", r1), re.findall("Zs:Z:(.*)\s?", r2))


def run():
    index = createGene()
    sample_index = index = createRead(index)
    index = kg_build_index.main(index)
    sample_index += hisatMap(index, sample_index)
    inspectMapping(sample_index)
    print(index, sample_index)


run()
