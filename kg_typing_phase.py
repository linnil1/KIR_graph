import os
from collections import defaultdict
from dataclasses import asdict

import pandas as pd
from pysam import VariantFile
from Bio import SeqIO

from namepipe import nt, NameTask
from kg_main import hisatKIRRessult, hisatMap, hisatTyping
from kg_utils import runDocker, threads
from kg_eval import EvaluateKIR
from kg_typping import readVariants, Variant, readAlleleLength
from kg_typping_linnil1 import AlleleTypping, readCNResult

import kg_build_index
from kg_msa_leftalign import fileLeftAlign

# if input_name.template_args[0] != "00":
#     return output_name


@nt
def createFastaIndex(index):
    # for GATK
    if os.path.exists(f"{index}_backbone.dict"):
        return index
    runDocker("gatk4", f""" \
        gatk CreateSequenceDictionary \
             -R {index}_backbone.fa
    """)
    return index


@nt
def msaLeftAlign(input_name):
    # kir_2100_2dl1s1.save.KIR2DL1S1.bam
    # -> kir_2100_2dl1s1.realign.save.KIR2DL1S1.bam
    output_temp = input_name.template.replace(".save", ".leftalign.save")
    output_name = input_name.replace(".save", ".leftalign.save")
    if os.path.exists(f"{output_name}.json"):
        return output_temp
    fileLeftAlign(input_name, output_name)
    return output_temp


@nt
def addRGBam(input_name):
    # add rg group in bam
    output_name = input_name + ".rg"
    if os.path.exists(f"{output_name}.bam"):
        return output_name

    runDocker("gatk4", f""" \
        gatk AddOrReplaceReadGroups \
             I={input_name}.bam \
             O={output_name}.bam \
             SORT_ORDER=coordinate \
             RGID={input_name} \
             RGLB=kir \
             RGPL=illumina \
             RGPU={input_name} \
             RGSM={input_name} \
             CREATE_INDEX=True
    """)
    return output_name


@nt
def sechuleCN(input_name):
    # split read by KIR gene

    # read CN
    # file_cn = "data4/linnil1_syn_30x_seed87.03.index_kir_2100_2dl1s1.mut01.hisatgenotype.errcorr.linnil1.cn_sam_depth_p75.tsv"
    file_cn = input_name.replace(".rg", ".hisatgenotype.errcorr.linnil1.cn_sam_depth_p75.tsv")
    gene_cn = readCNResult(file_cn)

    # per gene
    for gene_ref, cn in gene_cn.items():
        if not cn:
            continue
        output_name = input_name + "." + gene_ref.split("*")[0]
        if os.path.exists(f"{output_name}.bam"):
            continue

        # save copy number in txt
        open(f"{output_name}.txt", "w").write(str(cn))
        # spliit bam
        runDocker("gatk4", f""" \
            samtools view \
                     {input_name}.bam \
                     {gene_ref} \
                     -o {output_name}.bam \
        """)
        runDocker("gatk4", f""" \
            samtools index {output_name}.bam
        """)
    return input_name + ".{}"


@nt
def haplotypeCaller(input_name, index):
    # bam -> vcf
    output_name = input_name + ".hc"
    if not os.path.exists(f"{input_name}.txt"):
        return output_name
    if os.path.exists(f"{output_name}.g.vcf.gz"):
        return output_name

    cn = int(open(f"{input_name}.txt").read())
    runDocker("gatk4", f""" \
        gatk --java-options "-Xmx4g" HaplotypeCaller \
             -R {index}_backbone.fa \
             -I {input_name}.bam \
             -ploidy {cn} \
             -bamout {output_name}.bam \
             --create-output-bam-index true \
             -O {output_name}.g.vcf.gz
    """)
    return output_name


@nt
def whatHaps(input_name, old_input):
    # vcf 0/1 -> 0|1 (phased)
    output_name = input_name + ".wh"
    old_input_name = old_input.format(*input_name.template_args)

    if os.path.exists(f"{output_name}.g.vcf.gz"):
        return output_name

    cn = int(open(f"{old_input_name}.txt").read())
    runDocker("quay.io/biocontainers/whatshap:1.4--py37h96cfd12_0", f""" \
        whatshap polyphase \
        --threads {threads} \
        -B 1 \
        {input_name}.g.vcf.gz \
        {old_input_name}.bam \
        --indels \
        --ploidy {cn} \
        -o {output_name}.g.vcf.gz
    """)
    return output_name


@nt
def freebayes(input_name, index):
    # bam -> vcf
    output_name = input_name + ".fb"
    if not os.path.exists(f"{input_name}.txt"):
        return output_name
    if os.path.exists(f"{output_name}.g.vcf.gz"):
        return output_name

    cn = int(open(f"{input_name}.txt").read())
    runDocker("quay.io/biocontainers/freebayes:1.3.6--h346b5cb_1", f""" \
        freebayes \
        -f {index}_backbone.fa \
        --ploidy {cn} \
        --dont-left-align-indels \
        {input_name}.bam \
        --vcf {output_name}.g.vcf
    """)
    runDocker("samtools", f"bgzip {output_name}.g.vcf")
    return output_name


@nt
def atomizeAllele(input_name):
    # step 1:
    # normalize vcf allele (but not left-align)
    # step 2:
    # split strange variants int eh view of whatshap

    output_name1 = input_name + ".atomize"
    output_name2 = output_name1 + ".split"
    if os.path.exists(f"{output_name2}.g.vcf.gz"):
        return output_name2

    # step1
    runDocker("bcftools", f""" \
        bcftools norm -a \
            {input_name}.g.vcf.gz \
            -o {output_name1}.g.vcf.gz
    """)
    # split multi-allelic by '-m -both'

    # step2
    vcf = VariantFile(output_name1 + ".g.vcf.gz")
    vcf_out = VariantFile(output_name2 + ".g.vcf", 'w', header=vcf.header)

    for rec in vcf.fetch():
        if len(rec.ref) == 1:
            vcf_out.write(rec)
            continue

        r, a = rec.ref, rec.alts[0]
        if (a.startswith(r) or a.endswith(r) \
                or r.startswith(a) or r.endswith(a)):
            vcf_out.write(rec)
            continue

        assert len(a) == 1
        # case: 7066, 'TTTT', 'C'
        # -> T -> C
        i = 0
        new_rec = rec.copy()
        new_rec.pos += i
        new_rec.ref = r[i]
        new_rec.alts = [a[i]]
        vcf_out.write(new_rec)
        print(rec, "->", new_rec)

        # -> NTTT -> N
        new_rec = rec.copy()
        new_rec.pos += i
        new_rec.ref = "N" + r[1:]
        new_rec.alts = ["N"]
        vcf_out.write(new_rec)
        print(rec, "->", new_rec)

    vcf_out.close()
    runDocker("samtools", f"bgzip {output_name2}.g.vcf")
    return output_name2



def readVCF(file_vcf, variants_dict):
    # read vcf
    # yield (v, hap)
    # v: Variant, hap: it's belonging haplotype

    # init
    stat_all, stat_phase = 0, 0  # record level
    stat_alt_all, stat_alt_found = 0, 0  # alt level
    ploidy = None  # ploidy

    # read vcf records
    for rec in VariantFile(file_vcf).fetch():
        stat_all += 1
        sample = rec.samples.keys()[0]
        gt = rec.samples[sample]["GT"]
        if ploidy is None:
            ploidy = len(gt)
        else:
            assert ploidy == len(gt)

        # remove un-phase variants
        # but lower the accuracy (WHY)
        if not rec.samples[sample].phased:
            continue
        stat_phase += 1

        # format
        # chrom = KIR2DL1S1*BACKBONE 
        # pos = 2470
        # rec.ref = "AGT"
        # rec.alts = ('A', 'TGT')
        # gt = (1, 1, 2)
        # per alts
        for ai, alt in enumerate(rec.alts):
            if ai + 1 not in gt:
                continue
            stat_alt_all += 1
            v = None
            # SNP
            if len(alt) == 1 and len(rec.ref) == 1:
                v = Variant(pos=int(rec.pos) - 1, typ="single", val=alt, ref=rec.chrom)

            # after atomize + split -> this case can never occur
            elif len(alt) == len(rec.ref):
                assert False

            # deletion
            elif len(rec.ref) > len(alt):
                assert len(alt) == 1  # after atomize + split
                a_ref = rec.ref
                a_alt = alt
                k = 0

                # remove same left or right same base
                while len(a_alt) and a_ref[-1] == a_alt[-1]:
                    a_ref = a_ref[:-1]
                    a_alt = a_alt[:-1]
                while len(a_alt) and a_ref[0] == a_alt[0]:
                    a_ref = a_ref[1:]
                    a_alt = a_alt[1:]
                    k += 1

                # MSA is not left aligned: try each position -4 ~ 6
                for k in range(10):
                    v = Variant(pos=int(rec.pos) - 4 + k, typ="deletion", val=len(rec.ref) - len(alt), ref=rec.chrom)
                    if v in variants_dict:
                        break
                else:
                    print("Deletion not found", rec)
            else:
                # ignore insertion
                print(f"Insertion {rec}")

            # variant not found
            if v is None:
                continue
            if v not in variants_dict:
                print("\n\nNot found ", rec)
                print(asdict(v))
                continue

            # yield
            stat_alt_found += 1
            for gi, g in enumerate(gt):
                if g == ai + 1:
                    yield variants_dict[v], gi

    print(f"{stat_all=} {stat_phase}")
    print(f"{stat_alt_all=} {stat_alt_found}")
    # special case: vcf info
    yield None, ploidy


@nt
def alleleCalling(input_name, index):
    # vcf -> csv

    output_name = input_name + ".call"
    # if os.path.exists(f"{output_name}.csv"):
    #     return output_name
    # if str(input_name) != "data4/linnil1_syn_30x_seed87.09.data4_kir_2100_2dl1s1.leftalign.mut01.rg.KIR2DS3.fb.atomize.split.wh":
    #     return output_name

    # global variants
    variants = readVariants(index)
    variants_dict = {v: v for v in variants}

    # Gene
    # ref = input_name.split(".")[-3] + "*BACKBONE"
    ref = input_name.template_args[-1] + "*BACKBONE"
    print(ref)

    # read vcf
    haps_vars = defaultdict(list)
    for v, hap in readVCF(f"{input_name}.g.vcf.gz", variants_dict):
        if v is None:  # special case: vcf info
            for i in range(hap):
                if i not in haps_vars:
                    haps_vars[i] = []
            continue
        haps_vars[hap].append(v)

    if not len(haps_vars):
        haps_vars[0] = []
    # print(haps_vars)

    # calling
    called_alleles = []
    seq_len = readAlleleLength(index)

    for hap, variants in haps_vars.items():
        reads = []
        reads.extend([{
            'lp': [asdict(v)],
            'ln': [],
            'rp': [], 'rn': [],
        } for v in variants])
        if variants:
            assert variants[0].ref == ref
        reads.extend([{
            'lp': [],
            'ln': [asdict(v)],
            'rp': [], 'rn': [],
        } for v in set(i for i in variants_dict.keys() if i.ref == ref)- set(variants)])

        typing = AlleleTypping(reads, seq_len)
        result = typing.typingByNum(1)
        alleles = result['allele_name'][0]
        called_alleles.extend(alleles)  # must be length = 1

    print(called_alleles)
    open(f"{output_name}.csv", "w").write(",".join(called_alleles))
    return output_name


@nt
def mergeGeneCall(input_name):
    # genes' csv -> csv
    output_name = input_name.replace_wildcard("_merge_gene")
    # if os.path.exists(f"{output_name}.tsv"):
    #     return output_name

    called_alleles = []
    for name in input_name.get_input_names():
        called_alleles.extend(open(name + ".csv").read().split(','))

    df = pd.DataFrame([{
        'id': name.template_args[0],
        'name': input_name,
        'alleles': '_'.join(called_alleles),
    }])
    df.to_csv(f"{output_name}.tsv", index=False, sep="\t")
    return output_name


if __name__ == "__main__":
    answer = "linnil1_syn_30x_seed87"
    index = "data4/kir_2100_2dl1s1.save.{}"
    index = index >> msaLeftAlign >> NameTask(func=lambda i: i.template[:-3]) >> NameTask(func=kg_build_index.main)
    samples = "data4/linnil1_syn_30x_seed87.{}" >> hisatMap.set_args(index=str(index))

    # hisatgenotype filtering
    # samples = samples >> hisatTyping.set_args(index=str(index)) >> NameTask(func=lambda i: i + ".no_multi")

    # index >> createFastaIndex
    index = str(index)
    # >> haplotypeCaller.set_args(index=index) \
    # 
    samples = \
    samples >> addRGBam \
            >> sechuleCN
    samples >> freebayes.set_args(index=index) \
            >> atomizeAllele \
            >> whatHaps.set_args(str(samples)) \
            >> alleleCalling.set_args(index=index) \
            >> mergeGeneCall.set_depended(-1) \
            >> hisatKIRRessult.set_depended(-1).set_args(answer=answer)
