import pandas as pd
from namepipe import NamePath


def calcPossibleAlleleMean(name):
    dfs_possible = []
    for input_name in NamePath(name).get_input_names():
        df_possible = pd.read_csv(input_name + ".possible.tsv", sep="\t")
        df_possible["cn"] = df_possible.notna().sum(axis=1) - 3
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


def calcNovelAlleleMean(name):
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


if __name__ == "__main__":
    name = "data/linnil1_syn_s2022.{}.30x_s1031.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph"
           ".variant.noerrcorr.no_multi.depth.p75.CNgroup_b2_assume3DL3.pv_exonfirst_1.compare_sum.var_errcorr.top600"
    name = "data_real/hprc.{}.index_hs37d5.bwa.part_strict.index_kir_2100_withexon_ab_2dl1s1.leftalign.mut01.graph.trim"
           ".variant.noerrcorr.no_multi.depth.p75.CNgroup_b2_assume3DL3.pv_exonfirst_1.compare_sum.var_errcorr.top600"
    calcPossibleAlleleMean(name)
    calcNovelAlleleMean(name)
