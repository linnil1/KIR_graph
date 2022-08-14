"""
GRAPH-KIR: Main function
"""
import os
import argparse
from glob import glob
import pandas as pd

from gk_build_msa import buildKirMsa
from gk_build_index import msa2HisatReference, buildHisatIndex
from gk_hisat2 import hisatMap, extractVariantFromBam, readExons
from gk_cn import predictSamplesCN, loadCN
from gk_kir_typing import TypingWithPosNegAllele, TypingWithReport, Typing
from gk_plot import plotCN, plotReadMappingStat, showPlot, plotGeneDepths


def mergeAllele(allele_result_files: list[str], final_result_file: str):
    """ Merge allele calling result """
    df = pd.concat(pd.read_csv(f, sep="\t") for f in allele_result_files)
    df.to_csv(final_result_file, index=False, sep="\t")
    print(df)


def buildIndex(msa_type: str) -> tuple[str, str]:
    """
    Build graph data

    * MSA
    * HISAT2 formatted files (for hisat-build)
    * HISAT2 index

    Returns:
      HISAT2 files prefix
      HISAT2 index prefix
    """
    msa_index = f"index/kir_2100_{msa_type}"
    if not len(glob(msa_index + ".KIR*")):
        buildKirMsa(msa_type, msa_index)

    ref_index = msa_index + ".mut01"
    if not os.path.exists(ref_index + ".haplotype"):
        msa2HisatReference(msa_index, ref_index)

    index = ref_index + ".graph"
    if not os.path.exists(index + ".8.ht2"):
        buildHisatIndex(ref_index, index)
    return ref_index, index


def createParser() -> argparse.ArgumentParser:
    """ Setup arguments """
    parser = argparse.ArgumentParser(
            description='Run Graph KIR',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--r1", nargs="+",
                        default=["data/linnil1_syn_30x.00.read.1.fq"],
                        help="Read 1 fastq")
    parser.add_argument("--r2", nargs="+",
                        default=["data/linnil1_syn_30x.00.read.2.fq"],
                        help="Read 2 fastq")
    parser.add_argument("--msa-type",
                        default="ab_2dl1s1",
                        help="Type of MSA: merge, split, ab, ab_2dl1s1")
    parser.add_argument("--exon", action='store_true',
                        help="Select exon-only region of CN prediction")
    parser.add_argument("--cn-multi", action='store_true',
                        help="Do not remove multiple aligned reads for predicting CN")
    parser.add_argument("--cn-individually", action='store_true',
                        help="Predict CN separtely instead of cohort-wisely")
    parser.add_argument("--cn-cohort-model-output",  # require
                        default="data/linnil1_syn_30x.cohort.cn",
                        help="Cohort CN prediction model output path")
    parser.add_argument("--cn-select", default="p75",
                        help="Select 75 percentile(p75), mean(mean), median(median)"
                             " of depths")
    parser.add_argument("--cn-cluster", default="CNgroup",
                        help="CN prediction model (CNgroup, KDE)")
    parser.add_argument("--cn-provided", nargs="*",
                        help="Provided CN prediction in TSV")
    parser.add_argument("--allele-method", default="pv",
                        help="Max Likelihood by variants(pv) or report(report)")
    parser.add_argument("--result",  # require
                        default="data/linnil1_syn_30x.final.tsv",
                        help="Predicted alleles")
    return parser


def extractName(r1s: list[str], r2s: list[str]) -> list[tuple[str, str, str]]:
    """
    Extract the name for r1 and r2.

    The name is the longest common string between r1 and r2.

    Example:
      Input:
      * "data/linnil1_syn_30x.00.read.r1.fq" "data/linnil1_syn_30x.00.read.r2.fq"
      * "data/linnil1_syn_30x.01.read.r1.fq" "data/linnil1_syn_30x.01.read.r2.fq"

      Output's name:
      * "data/linnil1_syn_30x.00"
      * "data/linnil1_syn_30x.01"

    Returns:
      list of tuple[name, r1, r2]
    """
    data = []
    for r1, r2 in zip(r1s, r2s):
        name = ""
        for s1, s2 in zip(r1.split("."), r2.split(".")):
            if s1 == s2:
                if name:
                    name += "." + s1
                else:
                    name = s1
            else:
                break
        data.append((name, r1, r2))
    return data


def main(args: argparse.Namespace):
    """
    * Read Mapping
    * Extract variants
    * CN prediction
    * allele typing
    * merge
    """
    # Configuration
    data = extractName(args.r1, args.r2)

    # Main function
    cohort = not args.cn_individually
    ref_index, index = buildIndex(args.msa_type)

    exon_regions = {}
    if args.exon:
        exon_regions = readExons(ref_index)

    figs = []
    names = []
    bam_files = []
    cn_files = []
    suffix_cn_base = ".no_multi" if not args.cn_multi else ""

    for name, fq1, fq2 in data:
        suffix = "." + index.replace(".", "_").replace("/", "_")
        hisatMap(index, fq1, fq2, name + suffix + ".bam")
        name += suffix
        bam_files.append(name + ".bam")

        suffix = ".variant"
        extractVariantFromBam(ref_index, name + ".bam", name + suffix)
        name += suffix
        names.append(name)

        if not args.cn_provided and not cohort:
            suffix = f"{suffix_cn_base}.depth{'.exon' if args.exon else ''}"
            suffix_cn = f"{suffix}.{args.cn_select}.{args.cn_cluster}"
            predictSamplesCN([name + suffix_cn_base + ".bam"],
                             [name + suffix_cn + ".tsv"],
                             bam_selected_regions=exon_regions,
                             cluster_method=args.cn_cluster,
                             save_cn_model_path=name + suffix + ".json",
                             select_mode=args.cn_select)
            figs.extend(plotCN(name + suffix_cn + ".json"))
            figs.extend(plotGeneDepths(name + suffix + ".tsv"))
            cn_files.append(name + suffix_cn + ".tsv")
    figs.extend(plotReadMappingStat(bam_files, [fq1 for _, fq1, _ in data]))

    if not args.cn_provided and cohort:
        suffix = f"{suffix_cn_base}.depth{'.exon' if args.exon else ''}"
        suffix_cn = f"{suffix}.{args.cn_select}.cohort.{args.cn_cluster}"
        cohort_name = args.cn_cohort_model_output + suffix_cn
        predictSamplesCN([name + suffix_cn_base + ".bam" for name in names],
                         [name + suffix_cn + ".tsv" for name in names],
                         bam_selected_regions=exon_regions,
                         save_cn_model_path=cohort_name + ".json",
                         cluster_method=args.cn_cluster,
                         select_mode=args.cn_select)
        cn_files = [name + suffix_cn + ".tsv" for name in names]
        figs.append(plotCN(cohort_name + ".json"))

    if args.cn_provided:
        cn_files = args.cn_provided

    new_names = []
    for name, cn_file in zip(names, cn_files):
        result_name = name \
                      + ".cn_" + cn_file.replace("/", "_").replace(".", "_") \
                      + "." + args.allele_method

        if args.allele_method:
            t = TypingWithPosNegAllele(name + ".json")  # type: Typing
        else:
            t = TypingWithReport(name + ".em.json")
        cn = loadCN(cn_file)

        print(name)
        print(cn)
        called_alleles = t.typing(cn)
        print(called_alleles)

        t.save(result_name + ".json")
        df = pd.DataFrame([{
            'name': result_name,
            'alleles': "_".join(called_alleles),
        }])
        df.to_csv(result_name + ".tsv", sep="\t", index=False)
        new_names.append(result_name)

    mergeAllele([name + ".tsv" for name in new_names], args.result)
    return figs


def entrypoint():
    """ Command line entry function """
    parser = createParser()
    args = parser.parse_args()
    figs = main(args)
    showPlot(figs)


if __name__ == "__main__":
    entrypoint()
