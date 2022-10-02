"""
Raw depths -> gene depths -> gene's CN
"""
from itertools import chain
import pandas as pd

from .cn_model import CNgroup, KDEcut, Dist
from .utils import runDocker


def bam2Depth(file_bam: str, file_depth: str):
    """ Get read depth of all the position (via samtools depth) """
    runDocker("samtools", f"samtools depth -a {file_bam} -o {file_depth}")


def readSamtoolsDepth(depth_filename: str) -> pd.DataFrame:
    """ Read depths from samtools depths command (columns: gene, pos, depth) """
    df = pd.read_csv(depth_filename, sep="\t",
                     header=None, names=["gene", "pos", "depth"])
    return df


def selectSamtoolsDepth(df: pd.DataFrame,
                        ref_regions: dict[str, list[tuple[int, int]]]
                        ) -> pd.DataFrame:
    """ Depths outside regions are removed (Used for remove intron depth) """
    df_exon = []
    for gene, regions in ref_regions.items():
        for (start, end) in regions:
            df_exon.append(df[(df["gene"] == gene)
                           & (start <= df["pos"])
                           & (df["pos"] <= end)])
    return pd.concat(df_exon)


def aggrDepths(depths: pd.DataFrame, select_mode: str = "p75") -> pd.DataFrame:
    """ Depths of positions in each gene -> depths of gene """
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
              assume_3DL3_diploid: bool = False
              ) -> tuple[list[dict[str, int]], Dist]:
    """
    Depths of gene -> CN of gene

    The `assume_3DL3_diploid` assumption is useful when len(sample_gene_depths) == 1

    Parameters:
      sample_gene_depths: Gene's depth per sample
      assume_3DL3_diploid: Tuning the parameters by assuming KIR3DL3's CN is 2

    Returns:
      CN of gene per sample
      CN model
    """
    values = list(chain.from_iterable(map(lambda i: i['depth'], sample_gene_depths)))

    # cluster
    if cluster_method == "CNgroup":
        dist = CNgroup()
        dist.fit(values)
        if assume_3DL3_diploid:
            kir3dl3_depths = [
                float(gene_depths.loc[gene_depths['gene'] == "KIR3DL3*BACKBONE", "depth"])
                for gene_depths in sample_gene_depths]
            cn = dist.assignCN(kir3dl3_depths)
            if all(i == 1 for i in cn):
                dist.base /= 2

            cn = dist.assignCN(kir3dl3_depths)
            assert all(i == 2 for i in cn)

        # TODO: assume depth
        # if dist.base / assume_base > 1.7:
        #     dist.base /= 2

    elif cluster_method == "KDEcut":
        dist = KDEcut()  # type: ignore
        dist.fit(values)
    else:
        raise NotImplementedError

    sample_gene_cns = []
    for gene_depths in sample_gene_depths:
        sample_gene_cns.append(dict(zip(
            gene_depths['gene'],
            dist.assignCN(list(gene_depths['depth']))
        )))
    return sample_gene_cns, dist


def predictSamplesCN(samples_bam: list[str],
                     samples_cn: list[str],
                     bam_selected_regions: dict[str, list[tuple[int, int]]] = {},
                     save_cn_model_path: str | None = None,
                     assume_3DL3_diploid: bool = False,
                     select_mode: str = "p75",
                     cluster_method: str = "CNgroup"
                     ):
    """
    Read bamfile and predict CN per gene per sample

    Parameters:
      samples_bam: Bamfile path per sample
      samples_cn: The output CN file path per sample
      save_cn_model_path: Save the parameters of CN model into specific path
      bam_selected_regions:
        Use selected regions of read depths.
        Format: `dict[key=referce, value=list of tuple[start-position, end_position]]`
        Leave Empty to selected all regions
    """
    assert len(samples_bam) == len(samples_cn)

    # read bam -> depth per position -> depth per gene
    names = [bam.replace(".bam", "") for bam in samples_bam]
    sample_gene_depths = []
    for name in names:
        bam2Depth(name + ".bam", name + ".depth.tsv")
        depths = readSamtoolsDepth(name + ".depth.tsv")
        if bam_selected_regions:
            depths = selectSamtoolsDepth(depths, bam_selected_regions)
            depths.to_csv(name + ".depth.exon.tsv", index=False, sep="\t")
        gene_depths = aggrDepths(depths, select_mode=select_mode)
        sample_gene_depths.append(gene_depths)

    # depth per gene -> cn per gene
    cns, model = depthToCN(sample_gene_depths,
                           cluster_method=cluster_method,
                           assume_3DL3_diploid=assume_3DL3_diploid)

    if save_cn_model_path:
        model.save(save_cn_model_path)

    # cn per gene -> save
    for filename, cn in zip(samples_cn, cns):
        df = pd.DataFrame(list(cn.items()), columns=["gene", "cn"])
        df.to_csv(filename, index=False, sep="\t")


def loadCN(filename_cn: str) -> dict[str, int]:
    """
    Load CN data

    Return:
      key:  Reference name
      value:  CN on the reference
    """
    data = pd.read_csv(filename_cn, sep="\t", index_col=[0])
    return data.to_dict()['cn']
