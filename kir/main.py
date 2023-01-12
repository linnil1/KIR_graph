"""
Usage: kirpipe example_data/test.{}
"""
import argparse
from collections import defaultdict

import pandas as pd

from .kir_pipe import KirPipe
from .ping import PING
from .t1k import T1k
from .sakauekir import SakaueKir
from .kpi import KPI


def readArgument(factory: dict[str, KirPipe]) -> argparse.Namespace:
    """Read command line arguments"""
    parser = argparse.ArgumentParser(
        prog="KIR collections",
        description="Run all KIR tools",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "sample_name",
        help="samples name. Use {} to indicate wildcard.",
    )
    # parser.add_argument("--version", default="origin", help="IMGT-HLA version")
    parser.add_argument("--thread", default=4, help="Number of threads")
    parser.add_argument(
        "--tools",
        default=factory.keys(),
        nargs="+",
        help="KIR tools to execute",
    )
    parser.add_argument(
        "--ipd-version",
        default="2100",
        help="IPD-KIR database version (only works in some tools)"
    )
    parser.add_argument("--final-name", help="The name of final merged results")
    args = parser.parse_args()
    # print(args)
    return args


def groupAllele(result: pd.Series) -> pd.Series:  # type: ignore
    """A helper function of concatResult"""
    alleles = result["alleles"].split("_")
    gene = defaultdict(list)
    for allele in alleles:
        gene[allele.split("*")[0]].append(allele)
    return pd.Series(
        {
            **gene,
            "id": result["id"],
            "name": result["name"][:47],  # WARNING: Force to shorten the name
        }
    )


def concatResult(csv_files: list[str], output_name: str = "") -> pd.DataFrame:
    """Merge all KIR allele calling result"""
    df_list = [pd.read_csv(f + ".tsv", sep="\t", dtype=str) for f in csv_files]
    df = pd.concat(df_list)
    if output_name:
        df.to_csv(output_name, index=False, sep="\t")
    return df


def showResult(df: pd.DataFrame) -> pd.DataFrame:
    """Print result on screen"""
    new_df = df.apply(groupAllele, axis=1)
    new_df = new_df.melt(
        id_vars=["id", "name"],
        var_name="gene",
        value_name="allele",
    )
    with pd.option_context("display.max_rows", None):  # type: ignore
        print(new_df.groupby(["id", "gene", "name"])["allele"].apply(list))
    return new_df


def main() -> None:
    """Main entrypoint"""
    factory = {
        PING.name: PING(),
        PING.name + "-wgs": PING(version="wgs"),
        T1k.name: T1k(),
        SakaueKir.name: SakaueKir(),
        KPI.name: KPI(),
    }
    args = readArgument(factory)
    samples = args.sample_name

    results = []
    for tool in args.tools:
        module = factory[tool]
        module.setIPDVersion(args.ipd_version)
        module.setThreads(args.thread)
        result = module.runAll(samples)
        results.append(result)

    df = concatResult(results, output_name=args.final_name)
    showResult(df)


if __name__ == "__main__":
    main()
