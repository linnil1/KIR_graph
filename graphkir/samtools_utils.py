"""
Samtools utility functions for depth calculation
"""

import pandas as pd
from .external_tools import runTool


def bam2Depth(file_bam: str, file_depth: str, get_all: bool = True) -> None:
    """Get read depth of all the position (via samtools depth)"""
    if get_all:
        runTool("samtools", ["samtools", "depth", "-aa", file_bam, "-o", file_depth])
    else:
        runTool("samtools", ["samtools", "depth", file_bam, "-o", file_depth])


def readSamtoolsDepth(depth_filename: str) -> pd.DataFrame:
    """Read depths from samtools depths command (columns: gene, pos, depth)"""
    df = pd.read_csv(
        depth_filename, sep="\t", header=None, names=["gene", "pos", "depth"]
    )
    return df
