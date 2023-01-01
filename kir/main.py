"""
Usage: kirpipe example_data/test.{}
"""
import argparse
from typing import Type

from .kir_pipe import KirPipe
from .ping import PING
from .t1k import T1k
from .sakauekir import SakaueKir
from .kpi import KPI


def readArgument(factory: dict[str, Type[KirPipe]]) -> argparse.Namespace:
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
    args = parser.parse_args()
    # print(args)
    return args


def main() -> None:
    factory = {
        PING.name: PING,
        T1k.name: T1k,
        SakaueKir.name: SakaueKir,
        KPI.name: KPI,
    }
    args = readArgument(factory)
    samples = args.sample_name

    for tool in args.tools:
        module = factory[tool]
        result = module(threads=args.thread).runAll(samples)


if __name__ == "__main__":
    main()
