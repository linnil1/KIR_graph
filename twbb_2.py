import os
from glob import glob
from pathlib import Path
from functools import partial
import pandas as pd

from namepipe import nt, NameTask, compose, ConcurrentTaskExecutor
from graphkir.utils import runShell
from kg_main import buildMsaWithCds, leftAlignWrap, msa2HisatReferenceWrap, buildHisatIndexWrap, extractVariant, cnPredict, plotCNWrap


if __name__ == "__main__":
    data_folder = "data_twbb"
    index_folder = "index5"
    Path(data_folder).mkdir(exist_ok=True)
    Path(index_folder).mkdir(exist_ok=True)
    extract_exon = False

    # msa_index = index_folder >> buildKirMsaWrap.set_args("ab_2dl1s1") >> leftAlignWrap
    ref_index = compose([
        index_folder,
        partial(buildMsaWithCds, msa_type="ab_2dl1s1"),
        leftAlignWrap,
        msa2HisatReferenceWrap,
    ])
    index = ref_index >> buildHisatIndexWrap

    # mapping = "data_twbb/twbb.hg19.{}.part_merge.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim"
    # mapping = samples >> hisatMapWrap.set_args(index=str(index))
    mapping = "data_twbb/twbb.{}.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim"
    NameTask.default_executor = ConcurrentTaskExecutor()
    NameTask.default_executor.threads = 10
    variant = mapping >> nt(extractVariant).set_args(ref_index=str(ref_index))

    cn = variant >> nt(cnPredict).set_args(ref_index=str(ref_index), exon=extract_exon)  # .set_depended(0)

    # cn = variant >> cnPredict.set_args(ref_index=str(ref_index), exon=extract_exon).set_depended(0)

    # typing = variant >> kirTyping.set_args(cn, "pv_exonfirst") >> kirResult.set_args(answer=answer_folder).set_depended(0)
    # cn >> plotCNWrap.set_depended(0)
    runShell("stty echo opost")
