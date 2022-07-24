from pyhlamsa import msaio


def findDeletePos(seq):
    pos = 0
    length = 0
    for i in range(len(seq)):
        if seq[i] == "-":
            length += 1
        else:
            yield pos, length
            pos = i + 1
            length = 0
    if length:
        yield pos, length


def findShift(ref_seq, seq, pos, length):
    ref_seq_del = ref_seq[pos:pos+length]
    for i in range(length, 0, -1):
        if seq[pos - i:pos] == ref_seq_del[:i] == ref_seq_del[-i:]:
            return i
    return None


def applyShift(ref_seq, seq, pos, length, shift):
    # ACCACCACCACC
    # ACCTCC---ACC
    # ACC---TCCACC
    assert shift <= length
    new_seq = seq[:pos - shift]                        + seq[pos: pos + length] \
                                + seq[pos - shift:pos]                          + seq[pos + length:]
    assert new_seq.replace('-', '') == seq.replace('-', '')

    # print(f"{pos=} {shift=}")
    # print(ref_seq[pos - shift - 5: pos + length + 5])
    # print(    seq[pos - shift - 5: pos + length + 5])
    # print(new_seq[pos - shift - 5: pos + length + 5])
    return pos - shift, length, new_seq


def msaLeftAlign(msa):
    for name, seq in msa.alleles.items():
        ref_seq = msa.get_reference()[1]
        for pos, length in findDeletePos(seq):
            while True:
                shift = findShift(ref_seq, seq, pos, length)
                if shift is None:
                    break
                pos, length, seq = applyShift(ref_seq, seq, pos, length, shift)
                msa.alleles[name] = seq
    return msa


def fileLeftAlign(input_name, output_name):
    msa = msaio.load_msa(f"{input_name}.fa", f"{input_name}.json")
    msa = msaLeftAlign(msa)
    msaio.save_msa(msa, f"{output_name}.fa", f"{output_name}.json")


if __name__ == "__main__":
    gene = "KIR2DL1S1"
    msa = msaio.load_msa(f"data4/kir_2100_2dl1s1.save.{gene}.fa",
                         f"data4/kir_2100_2dl1s1.save.{gene}.json")
