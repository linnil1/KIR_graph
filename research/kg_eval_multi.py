from typing import Iterable

import pandas as pd
from namepipe import NamePath

from kg_eval import readAnswerAllele, compareSample, MatchType, compareCohort


def calcPossibleAlleleMean(name: str) -> None:
    dfs_possible = []
    for input_name in NamePath(name).get_input_names():
        df_possible = pd.read_csv(input_name + ".possible.tsv", sep="\t")
        df_possible["cn"] = map(len, map(extractAlleleFromPossibleFormat(df_possible)))
        # print(df_possible)
        dfs_possible.append(
            df_possible.groupby(["gene"], as_index=False).agg(
                {"rank": "size", "cn": "mean"}
            )
        )

    df = pd.concat(dfs_possible)
    print(df.groupby(["gene", "cn"])["rank"].mean())
    print(df.groupby(["gene"])["rank"].mean())
    print(df.groupby(["cn"])["rank"].mean())


def calcNovelAlleleMean(name: str) -> None:
    dfs_variant = []
    for input_name in NamePath(name).get_input_names():
        df_variant = pd.read_csv(input_name + ".novel.variant.tsv", sep="\t")
        df_variant["sample"] = input_name.template_args[-1]
        dfs_variant.append(df_variant)
        # df_variant = df_variant[df_variant["skip"] == False]
        # print(df_variant)
        # len(df_variant['skip'])
    df_variant = pd.concat(dfs_variant)
    df_variant["not_skip"] = ~df_variant["skip"]
    df_variant["allele_id"] = df_variant.groupby(["sample", "allele"]).ngroup()

    df_allele_variant_num = df_variant.groupby(["allele_id"]).agg(
        {"not_skip": "sum", "gene": "first", "sample": "first"}
    )
    df_variant_mean = pd.concat(
        [
            df_allele_variant_num.groupby("gene").agg({"not_skip": "mean"}),
            df_allele_variant_num.groupby(["sample", "gene"])
            .agg({"not_skip": "sum"})
            .groupby("gene")
            .agg({"not_skip": "mean"}),
        ],
        axis=1,
    )
    df_variant_mean.set_axis(
        ["Mean_per_allele", "Mean_per_sample"], axis="columns", inplace=True
    )
    print(df_variant_mean)


def extractAlleleFromPossibleFormat(df: pd.DataFrame) -> Iterable[list[str]]:
    allele_num = set(map(str, range(10)))
    for _, txt in df.iterrows():
        yield [v for k, v in dict(txt).items() if type(v) is str and k in allele_num]


def findBestSet(
    answer_alleles: list[str], allele_sets: Iterable[list[str]]
) -> list[str]:
    best_set = []
    best_score = 0
    for allele_set in allele_sets:
        data = compareSample(answer_alleles, allele_set)
        score = 0
        for i in data:
            if i.match_type == MatchType.MATCH7:
                score += 7
            if i.match_type == MatchType.MATCH5:
                score += 4
            if i.match_type == MatchType.MATCH3:
                score += 3
            if i.match_type == MatchType.MATCHGENE:
                score += 1
        # print(allele_set, score)
        if score > best_score:
            best_score = score
            best_set = allele_set
    return best_set


def evalMultiAllele(sample: str, answer: str) -> None:
    """Find the best allele set in possible sets and compare to answer"""
    output_name = NamePath(sample).replace_wildcard("_merge_answer_best")
    # if Path(output_name + ".tsv").exists():
    cohort_answer = readAnswerAllele(answer)
    cohort_best_sample = {}
    for name in NamePath(sample).list_names():
        sid = name.template_args[0]
        df = pd.read_csv(name + ".tsv", sep="\t")
        alleles = []
        for gene in sorted(set(df["gene"])):
            df_gene = df[df["gene"] == gene]
            allele_sets = extractAlleleFromPossibleFormat(df_gene)
            allele_set = findBestSet(cohort_answer[sid], allele_sets)
            alleles.extend(allele_set)
        cohort_best_sample[sid] = alleles

    compareCohort(
        cohort_answer,
        cohort_best_sample,
    )

    df = pd.DataFrame(
        [
            {"id": sid, "name": sample.format(sid), "alleles": "_".join(alleles)}
            for sid, alleles in cohort_best_sample.items()
        ]
    )
    df.to_csv(output_name + ".tsv", sep="\t", index=False)
    print(output_name)


if __name__ == "__main__":
    name = (
        "data/linnil1_syn_s2022.{}.30x_s1031.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph"
        + ".variant.noerrcorr.no_multi.depth.p75.CNgroup_b2_assume3DL3.pv_exonfirst_1.compare_sum.var_errcorr.top600"
    )
    name = (
        "data_real/hprc.{}.index_hs37d5.bwa.part_strict.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim"
        + ".variant.noerrcorr.no_multi.depth.p75.CNgroup_b2_assume3DL3.pv_exonfirst_1.compare_sum.var_errcorr.top600"
    )
    calcPossibleAlleleMean(name)
    calcNovelAlleleMean(name)
    exit()
    name = (
        "data/linnil1_syn_s2022.{}.30x_s1031.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph"
        + ".variant.noerrcorr.no_multi.depth.p75.CNgroup_b2_assume3DL3.pv.compare_sum.var_errcorr.top600.possible"
    )
    name = (
        "data/linnil1_syn_s2022.{}.30x_s1031.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph"
        + ".variant.noerrcorr.no_multi.depth.p75.CNgroup_b2_assume3DL3.pv_exonfirst_1.compare_sum.var_errcorr.top600.possible"
    )
    answer = "linnil1_syn/linnil1_syn_s2022_summary.csv"
    evalMultiAllele(name, answer)
