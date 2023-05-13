"""
Raw depths -> gene depths -> gene's CN
"""
from typing import Any
from itertools import chain
import json

import pandas as pd

from .utils import runDocker, NumpyEncoder, logger
from .cn_model import CNgroup, KDEcut, Dist


def bam2Depth(file_bam: str, file_depth: str, get_all: bool = True) -> None:
    """Get read depth of all the position (via samtools depth)"""
    if get_all:
        runDocker("samtools", f"samtools depth -aa {file_bam} -o {file_depth}")
    else:
        runDocker("samtools", f"samtools depth     {file_bam} -o {file_depth}")


def readSamtoolsDepth(depth_filename: str) -> pd.DataFrame:
    """Read depths from samtools depths command (columns: gene, pos, depth)"""
    df = pd.read_csv(
        depth_filename, sep="\t", header=None, names=["gene", "pos", "depth"]
    )
    return df


def selectSamtoolsDepth(
    df: pd.DataFrame, ref_regions: dict[str, list[tuple[int, int]]]
) -> pd.DataFrame:
    """Depths outside regions are removed (Used for remove intron depth)"""
    df_exon = []
    for gene, regions in ref_regions.items():
        for start, end in regions:
            df_exon.append(
                df[(df["gene"] == gene) & (start <= df["pos"]) & (df["pos"] <= end)]
            )
    return pd.concat(df_exon)


def aggrDepths(depths: pd.DataFrame, select_mode: str = "p75") -> pd.DataFrame:
    """Depths of positions in each gene -> depths of gene"""
    if select_mode == "median":
        gene_depths = depths.groupby(by="gene", as_index=False)["depth"].median()
    elif select_mode == "mean":
        gene_depths = depths.groupby(by="gene", as_index=False)["depth"].mean()
    elif select_mode == "p75":
        gene_depths = depths.groupby(by="gene", as_index=False)["depth"].quantile(0.75)
    else:
        raise NotImplementedError
    return gene_depths  # type: ignore


def depthToCN(
    sample_gene_depths: list[dict[str, float]],
    cluster_method: str = "CNgroup",
    cluster_method_kwargs: dict[str, Any] = {},
    assume_3DL3_diploid: bool = False,
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
    values = list(chain.from_iterable(map(lambda i: i.values(), sample_gene_depths)))
    logger.info(f"[CN] Predict copy number by {cluster_method} with data size {len(values)}")

    # cluster
    if cluster_method == "CNgroup":
        dist = CNgroup()
        if cluster_method_kwargs:
            dist = CNgroup.setParams(dist.getParams() | cluster_method_kwargs)
        logger.debug(f"[CN] Parameters before: {dist.getParams()}")
        dist.fit(values)
        if assume_3DL3_diploid:
            kir3dl3_depths = [
                float(gene_depths["KIR3DL3*BACKBONE"])
                for gene_depths in sample_gene_depths
            ]
            cn = dist.assignCN(kir3dl3_depths)
            if all(i == 1 for i in cn):
                logger.debug("[CN] Assume 3DL3 cn=2")
                assert isinstance(dist.base, float)
                dist.base *= 1 / 2
            # if all(i == 3 for i in cn):
            #     logger.debug("[CN] Assume 3DL3 cn=2")
            #     assert isinstance(dist.base, float)
            #     dist.base *= 3 / 2
            if all(i == 4 for i in cn):
                logger.debug("[CN] Assume 3DL3 cn=2")
                assert isinstance(dist.base, float)
                dist.base *= 4 / 2

            cn = dist.assignCN(kir3dl3_depths)
            assert all(i == 2 for i in cn)
        logger.info(f"[CN] {cluster_method} base = {dist.base}")
        # logger.debug(f"[CN] Parameters after:  {dist.getParams()}")

        # TODO: assume depth
        # if dist.base / assume_base > 1.7:
        #     dist.base /= 2

    elif cluster_method == "kde":
        dist = KDEcut()  # type: ignore
        dist.fit(values)
        logger.info(f"[CN] {cluster_method} cut = {dist.local_min}")  # type: ignore
        # logger.debug(f"[CN] Parameters after:  {dist.getParams()}")
    else:
        raise NotImplementedError

    sample_gene_cns = []
    for gene_depths in sample_gene_depths:
        genes, depths = zip(*gene_depths.items())
        sample_gene_cns.append(dict(zip(genes, dist.assignCN(depths))))  # type: ignore
    return sample_gene_cns, dist


def filterDepth(
    depth_file: str,
    filtered_depth_file: str,
    bam_selected_regions: dict[str, list[tuple[int, int]]] = {},
) -> None:
    """
    Bam -> tsv of read depth

    Parameters:
      bam_selected_regions:
        Use selected regions of read depths.
        Format: `dict[key=referce, value=list of tuple[start-position, end_position]]`
        Leave Empty to selected all regions
    """
    logger.debug("[Graph] exon: {bam_selected_regions}")
    depths = readSamtoolsDepth(depth_file)
    depths = selectSamtoolsDepth(depths, bam_selected_regions)
    depths.to_csv(filtered_depth_file, header=False, index=False, sep="\t")


def predictSamplesCN(
    samples_depth_tsv: list[str],
    samples_cn: list[str],
    save_cn_model_path: str | None = None,
    assume_3DL3_diploid: bool = False,
    select_mode: str = "p75",
    per_gene: bool = False,
    cluster_method: str = "CNgroup",
    cluster_method_kwargs: dict[str, Any] = {},
) -> None:
    """
    Read depth tsv and predict CN per gene per sample

    Parameters:
      samples_depth_tsv: Depths per sample (in tsv)
      samples_cn: The output CN file path per sample
      save_cn_model_path: Save the parameters of CN model into specific path
    """
    assert len(samples_depth_tsv) == len(samples_cn)

    # read bam -> depth per position -> depth per gene
    sample_gene_depths = []
    for depth_file in samples_depth_tsv:
        logger.info(f"[CN] Select {select_mode} of depths per gene ({depth_file})")
        df = aggrDepths(readSamtoolsDepth(depth_file), select_mode=select_mode)
        df["depth_file"] = depth_file
        sample_gene_depths.append(df)

    # TODO: If normalized needed, write here.
    logger.info(f"[CN] Predict CN from {len(sample_gene_depths)} samples")
    depths_dict = [dict(zip(i["gene"], i["depth"])) for i in sample_gene_depths]
    if not per_gene:
        # depth per gene -> cn per gene
        cns, model = depthToCN(
            depths_dict,
            cluster_method=cluster_method,
            cluster_method_kwargs=cluster_method_kwargs,
            assume_3DL3_diploid=assume_3DL3_diploid,
        )
        model.raw_df = [df.to_dict() for df in sample_gene_depths]
        if save_cn_model_path:
            model.save(save_cn_model_path)
    else:
        # concat sample with id
        file_index = {name: i for i, name in enumerate(samples_depth_tsv)}
        df_depths = pd.concat(sample_gene_depths)
        df_depths["gene_sampleid"] = df_depths["gene"] + "-" + df_depths["depth_file"]

        # save per gene cn
        cns = [{} for _ in range(len(sample_gene_depths))]
        cns_model = []
        for gene in sorted(set(df_depths["gene"])):
            logger.info(f"[CN] Predict per gene: {gene}")
            # extract same gene and use the fake name
            # bcz depthToCN use name as key, it'll overwrite
            gene_depths = df_depths[df_depths["gene"] == gene]
            gene_cns, gene_model = depthToCN(
                [dict(zip(gene_depths["gene_sampleid"], gene_depths["depth"]))],
                cluster_method=cluster_method,
                cluster_method_kwargs=cluster_method_kwargs,
                # assume_3DL3_diploid=assume_3DL3_diploid,
            )
            # print(gene_cns)
            # save back to per sample
            gene_model.raw_df = [gene_depths.to_dict()]
            cns_model.append((gene, gene_model))
            for gene_and_id in gene_cns[0]:
                i = file_index[gene_and_id.split("-")[1]]
                cns[i][gene] = gene_cns[0][gene_and_id]

        if save_cn_model_path:
            data = []
            for gene, model in cns_model:
                data.append(model.getParams())
                data[-1]["gene"] = gene
                json.dump(data[-1], open(save_cn_model_path + f".{gene}.json", "w"), cls=NumpyEncoder)
            json.dump(data, open(save_cn_model_path, "w"), cls=NumpyEncoder)

    # cn per gene -> save
    for filename, cn, depths in zip(samples_cn, cns, depths_dict):
        df1 = pd.DataFrame(list(cn.items()), columns=["gene", "cn"])
        df2 = pd.DataFrame(list(depths.items()), columns=["gene", "depth"])
        df = df1.merge(df2, on="gene")
        df.to_csv(filename, index=False, sep="\t")


def loadCN(filename_cn: str) -> dict[str, int]:
    """
    Load CN data

    Return:
      key:  Reference name
      value:  CN on the reference
    """
    data = pd.read_csv(filename_cn, sep="\t", index_col=[0])
    return dict(data.to_dict()["cn"])
