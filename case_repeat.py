from pyHLAMSA import Genemsa
from pprint import pprint
import os
from Bio import AlignIO


index = "index/kir_2100_merge"
gene = "KIR"
# index = "index/kir_2100_raw"
# gene = "KIR3DL2"


def run_dk(image, cmd):
    """ run docker container """
    run("docker run -it --rm --security-opt label=disable -u root -w /app "
        f"-v $PWD:/app {image} {cmd}")


def run(cmd):
    """ wrap os.system """
    print(cmd)
    os.system(cmd)


def seqSplit(seq, sep):
    seqs = []
    for s in seq:
        s_split = s.split(sep)
        s_split[1:] = [sep + i for i in s_split[1:]]
        seqs.extend(s_split)
    return seqs


def writeSeqs(seqs, filename):
    with open(filename, "w") as f:
        for seq in seqs:
            f.write(">" + seq['id'] + "\n")
            f.write(seq['seq'] + "\n")


def addPosition(seq_list):
    pos = 0
    seqs = []
    for seq in seq_list:
        seqs.append({
            'id': f"{gene}*BB-{pos:05d}",
            'seq': seq,
        })
        pos += len(seq)
    return seqs


def muscle(filename):
    run_dk("quay.io/biocontainers/muscle:5.1--h9f5acd7_1",
           f"muscle -align {filename}.fa -threads 2"
           f"      -output {filename}.muscle.fa")
    return AlignIO.read(filename+ ".muscle.fa", "fasta")


def selectUniqueName(names, sep="*"):
    s = set()
    for i in sorted(names):
        if sep == "*":
            key = i.split("*")[0]
        elif type(sep) is int:
            key = i.split("*")[0] + "*" + i.split("*")[1][:sep]
        else:
            key = i
        if key not in s:
            s.add(key)
            yield i
 


if __name__ == "__main__":
    msa = Genemsa.load_msa(f"{index}.{gene}.save.fa",
                           f"{index}.{gene}.save.json")
    # get repeat seqs
    seq_list = [msa.get(f"{gene}*BACKBONE")[:5000]]
    seq_list = seqSplit(seq_list, "AATATGG")
    seq_list = seqSplit(seq_list, "ATAGGG")
    seq_list = seqSplit(seq_list, "GATA")
    seq_list = seqSplit(seq_list, "GATC")
    seq_list = seqSplit(seq_list, "GAGAGGAA")
    seq_list = seqSplit(seq_list, "GAGTGAG")
    seq_list = seqSplit(seq_list, "GATGTG")
    seq_list = seqSplit(seq_list, "GGTAG")
    seq_list = seqSplit(seq_list, "GGGAGCG")
    seq_list = seqSplit(seq_list, "GGGAGTGCG")
    seq_list = seqSplit(seq_list, "GGGGAGA")
    seq_list = seqSplit(seq_list, "GTAAT")
    seq_list = seqSplit(seq_list, "GTGAT")
    seq_list = seqSplit(seq_list, "GTTAT")
    seq_list = addPosition(seq_list)
    # for i in seq_list:
    #     print(i['seq'])

    # reserve similar length
    seq_list = [i for i in seq_list if 18 <= len(i['seq']) <= 24]
    pprint(seq_list)

    # run muscle
    file_repeat = f"repeaet.{gene}"
    writeSeqs(seq_list, file_repeat + ".fa")
    msa_align = muscle(file_repeat)

    # print aligned repeat
    msa_align = Genemsa.from_MultipleSeqAlignment(msa_align)
    msa_align.append(f"{gene}*BB-con", msa_align.get_consensus(include_gap=True))
    msa_align.set_reference(f"{gene}*BB-con")
    print(msa_align.sort_name().format_alignment_diff())

    # print region of repeat of all alleles
    # alleles_show = list(selectUniqueName(msa.get_sequence_names(), sep=3))[::4]
    alleles_show = list(selectUniqueName(msa.get_sequence_names(), sep=5))[::4]
    msa = msa.select_allele(alleles_show)
    for name in sorted(msa_align.get_sequence_names()):
        if "-con" in name:
            continue
        length = len(msa_align.get(name).replace("-", ""))
        pos = int(name.split('-')[-1])
        print(msa[pos:pos+length].format_alignment_diff())
