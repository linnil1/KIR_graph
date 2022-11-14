"""
Raw depths -> gene depths -> gene's CN
"""
import json
from itertools import chain
import pandas as pd

from .cn_model import CNgroup, KDEcut, Dist
from .utils import runDocker, NumpyEncoder


def bam2Depth(file_bam: str, file_depth: str, get_all=True):
    """ Get read depth of all the position (via samtools depth) """
    if get_all:
        runDocker("samtools", f"samtools depth -a {file_bam} -o {file_depth}")
    else:
        runDocker("samtools", f"samtools depth    {file_bam} -o {file_depth}")


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


def filterDepth(depth_file: str,
                filtered_depth_file: str,
                bam_selected_regions: dict[str, list[tuple[int, int]]] = {}):
    """
    Bam -> tsv of read depth

    Parameters:
      bam_selected_regions:
        Use selected regions of read depths.
        Format: `dict[key=referce, value=list of tuple[start-position, end_position]]`
        Leave Empty to selected all regions
    """
    depths = readSamtoolsDepth(depth_file)
    depths = selectSamtoolsDepth(depths, bam_selected_regions)
    depths.to_csv(filtered_depth_file, header=False, index=False, sep="\t")


def predictSamplesCN(samples_depth_tsv: list[str],
                     samples_cn: list[str],
                     save_cn_model_path: str | None = None,
                     assume_3DL3_diploid: bool = False,
                     select_mode: str = "p75",
                     per_gene: bool = False,
                     cluster_method: str = "CNgroup",
                     ):
    """
    Read depth tsv and predict CN per gene per sample

    Parameters:
      samples_depth_tsv: Depths per sample (in tsv)
      samples_cn: The output CN file path per sample
      save_cn_model_path: Save the parameters of CN model into specific path
    """
    assert len(samples_depth_tsv) == len(samples_cn)

    # read bam -> depth per position -> depth per gene
    sample_gene_depths = [
        aggrDepths(readSamtoolsDepth(depth_file), select_mode=select_mode)
        for depth_file in samples_depth_tsv
    ]
    # TODO: If normalized needed, write here.

    if not per_gene:
        # depth per gene -> cn per gene
        cns, model = depthToCN(sample_gene_depths,
                               cluster_method=cluster_method,
                               assume_3DL3_diploid=assume_3DL3_diploid)

        if save_cn_model_path:
            model.save(save_cn_model_path)
    else:
        # concat sample with id
        for i, df in enumerate(sample_gene_depths):
            df['gene_sampleid'] = df['gene'] + "-" + str(i)
        depths = pd.concat(sample_gene_depths)

        # save per gene cn
        cns = [{} for _ in range(len(sample_gene_depths))]
        cns_model = []
        for gene in sorted(set(depths["gene"])):
            # print(gene)
            # extract same gene and use the fake name
            # bcz depthToCN use name as key, it'll overwrite
            gene_depths = depths[depths["gene"] == gene]
            gene_depths["gene"] = gene_depths["gene_sampleid"]
            gene_cns, gene_model = depthToCN([gene_depths],
                                             cluster_method=cluster_method,
                                             assume_3DL3_diploid=assume_3DL3_diploid)
            # print(gene_cns)
            # save back to per sample
            cns_model.append((gene, gene_model))
            for gene_and_id in gene_cns[0]:
                i = int(gene_and_id.split("-")[1])
                cns[i][gene] = gene_cns[0][gene_and_id]

        if save_cn_model_path:
            data = []
            for gene, model in cns_model:
                data.append(model.getParams())
                data[-1]["gene"] = gene
            json.dump(data, open(save_cn_model_path, "w"), cls=NumpyEncoder)

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
