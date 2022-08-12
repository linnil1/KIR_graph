"""
Raw depths -> gene depths -> gene's CN
"""
import pandas as pd
from gk_cn_model import CNgroup, KDEcut, Dist
from gk_utils import runDocker


def readSamtoolsDepth(depth_filename: str) -> pd.DataFrame:
    """ Read depths from samtools depths command """
    df = pd.read_csv(depth_filename, sep="\t",
                     header=None, names=["gene", "pos", "depth"])
    return df


def selectSamtoolsDepthInExon(df: pd.DataFrame,
                              exon_region: dict[str, list[tuple[int, int]]]
                              ) -> pd.DataFrame:
    """ intron depths are removed """
    df_exon = []
    for gene, regions in exon_region.items():
        for (start, end) in regions:
            df_exon.append(df[(df["gene"] == gene)
                           & (start <= df["pos"])
                           & (df["pos"] <= end)])
    return pd.concat(df_exon)


def aggrDepths(depths: pd.DataFrame, select_mode: str = "p75") -> pd.DataFrame:
    """ depths of positions -> depths of gene -> CN of gene """
    # choose depths
    if select_mode == "median":
        gene_depths = depths.groupby(by="gene", as_index=False)['depth'].median()
    elif select_mode == "mean":
        gene_depths = depths.groupby(by="gene", as_index=False)['depth'].mean()
    elif select_mode == "p75":
        gene_depths = depths.groupby(by="gene", as_index=False)['depth'].quantile(.75)

    else:
        raise NotImplementedError
    return gene_depths  # type: ignore


def depthToCN(sample_gene_depths: list[pd.DataFrame],
              cluster_method: str = "CNgroup",
              assume_3DL3_diploid=False) -> list[dict[str, int]]:
    """
    Depths of gene -> CN of gene
    The assumption is usful when len(sample_gene_depths) == 1
    """
    values = list(chain.from_iterable(map(lambda i: i['depth'], sample_gene_depths)))

    # cluster
    if cluster_method == "CNgroup":
        dist = CNgroup()  # type: Dist
        # assume_base=self.assume_depth)
        dist.fit(values)
        if assume_3DL3_diploid:
            kir3dl3_id = gene_depths['gene'] == "KIR3DL3"
            kir3dl3_depths = [float(gene_depths.loc[kir3dl3_id, "depth"])
                              for gene_depths in sample_gene_depths]
            cn = dist.assignCN(kir3dl3_depths)
            if all(i == 1 for i in cn):
                dist.base /= 2

            cn = dist.assignCN(kir3dl3_depths)
            assert all(i == 2 for i in cn)
        # cg.assume_3DL3_diploid = True
        # gene_cn = cg.predictKIRCN(gene_depth_dict)

    elif cluster_method == "KDEcut":
        dist = KDEcut()
        dist.fit(values)
        # gene_cn = dist.predictCN(gene_depth_dict)
    else:
        raise NotImplementedError

    # TODO
    # self.figs.extend(dist.plot())
    sample_gene_cns = []
    for gene_depths in sample_gene_depths:
        sample_gene_cns.append(dict(zip(
            gene_depths['gene'],
            dist.assignCN(list(gene_depths['depth']))
        )))
    return sample_gene_cns


def bam2Depth(file_bam: str, file_depth: str):
    """ Get read depth of all the position (via samtools depth) """
    runDocker("samtools", f"samtools depth -a {file_bam} -o {file_depth}")


def predictSamplesCN(samples_bam: list[str]) -> list[str]:
    """ """
    names = [bam.replace(".bam", "") for bam in samples_bam]

    sample_gene_depths = []
    for name in names:
        # bam2Depth(name + ".bam", name + ".depth.tsv")
        depths = readSamtoolsDepth(name + ".depth.tsv")
        gene_depths = aggrDepths(depths)
        sample_gene_depths.append(gene_depths)

    cns = depthToCN(sample_gene_depths)
    output_names = []
    for name, cn in zip(names, cns):
        print(name, cn)
        df = pd.DataFrame(list(cn.items()), columns=["gene", "cn"])
        print(name, df)
        df.to_csv(name + ".depth.cn_p75.tsv", index=False, sep="\t")
        output_names.append(name + ".depth.cn_p75.tsv")

    return output_names


if __name__ == "__main__":
    """
    name = "data/linnil1_syn_30x.00.index_kir_2100_ab_2dl1s1_muscle_mut01_graph.variant"
    bam2Depth(name + ".no_multi.bam",
              name + ".no_multi.depth.tsv")
    depths = readSamtoolsDepth(name + ".no_multi.depth.tsv")
    gene_depths = aggrDepths(depths)
    print(gene_depths)
    cn = depthToCN([gene_depths])
    print(cn[0])

    df = pd.DataFrame(cn[0].items(), columns=["gene", "cn"])
    print(df)
    df.to_csv(name + ".variant.no_multi.depth.cn_p75.tsv", index=False, sep="\t")
    """
    predictSamplesCN([name + ".no_multi.bam"])
