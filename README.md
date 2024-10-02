# Graph-KIR

Graph-KIR is a tool for KIR (Killer Immunoglobulin-like Receptor) typing using short read FASTQ files.

Paper Link: (Not yet published)

Biorxiv: https://doi.org/10.1101/2023.11.29.568665

github tag: v1.0


This repo contains two main programs:

1. `graphkir` - Main Typing Tool
  
    `graphkir` reads FASTQ files, both from CSV or directly via command-line arguments.
    It outputs copy number estimations in a CSV file called `cohort.cn.tsv` and 
    allele typing results in `cohort.allele.tsv` by default.
    More details about its algorithm and concept can be found in the paper.

2. `kirpipe` - KIR Typing Pipeline

    `kirpipe` is an aggregation tool that automates the KIR typing pipeline.
    It includes five published tools: `graphkir`, `PING`, `Sakaue's KIR`, `T1K`, and `KIR*KPI`.

    (Note: Currently, `kirpipe` requires podman or docker to execute)


## Requirements

Before using Graph-KIR, ensure you meet these requirements:

* Python >= 3.10

You have the option to use one of the following containerization tools:

* podman 
* docker
* singularity

You can choose to use it by specifying the `engine`. i.e. `--engine podman`.

If none of these containerization tools are installed, you can run Graph-KIR locally `--engine local`.
However, you'll need to install the following external packages:

* MUSCLE >= 5.1 (required only for index building stage)
* HISAT2 >= 2.2.1
* samtools >= 1.15.1
* BWA-MEM >= 0.7.17 (needed only for the WGS extraction stage)
* wget (necessary for downloading hs37d5 in the WGS extraction stage)


## Usage (Main)

Download the pre-built Graph-KIR index:
``` bash
wget https://graphkir.c4lab.tw/download/example_index.tar.gz
tar xvf example_index.tar.gz
# If kirpipe is used, rename it
# ln -s example_index graphkir_alpha
```

Install Graph-KIR:
``` bash
git clone https://github.com/linnil1/KIR_graph
cd KIR_graph
pip install .
graphkir --help
```

Run Graph-KIR (If the index does not exist, it will be auto-built):

``` bash
graphkir \
    --thread 2 \
    --r1 example/test00.read1.fq.gz \
    --r2 example/test00.read2.fq.gz \
    --r1 example/test01.read1.fq.gz \
    --r2 example/test01.read2.fq.gz \
    --index-folder example_index \
    --output-folder example_data \
    --output-cohort-name example_data/cohort
```

Or, if you have an input CSV file (e.g., `cohort.csv`) containing the list of samples:

``` bash
graphkir \
    --thread 2 \
    --input-csv example/cohort.csv \
    --index-folder example_index \
    --allele-method exonfirst \
    --output-cohort-name example_data/cohort \
    --log-level DEBUG
```

The CSV should have four columns:

* `name`: The output prefix of the sample.
* `r1` and `r2`: Paths to the fastq files.
* `cnfile`: You can assign a copy number file for the sample. Leave it empty for Graph-KIR to assign automatically.

``` csv
name,r1,r2,cnfile
example_data/linnil1.00,example/test00.read1.fq.gz,example/test00.read2.fq.gz,example/test00.assigned.cn.tsv
example_data/linnil1.01,example/test01.read1.fq.gz,example/test01.read2.fq.gz,
```

The final result that includes all the samples are aggrate into one file with prefix `output-cohort-name`.
In the above sample, `example_data/cohort.cn.tsv` and `example_data/cohort.allele.tsv` are generated.

Some useful arguments include:
* `--cn-cohort`: Estimate copy number while considering the entire cohort.
* `--cn-3dl3-not-diploid`: Do not assume that the copy number of 3DL3 is equal to 2.
* `--allele-strategy exonfirst`: Perform typing using the exon part of reads instead of the entire sequence.
* You can manually assess the copy number estimation results using the `--plot` option.
* Adjust the distribution deviation with the `--cn-dist-dev` argument, for example, `--cn-dist-dev 0.06`.


## Usage (`kirpipe` pipeline for other KIR tools)
```
ln -s ../example/test00.read1.fq.gz example_data/test.00.read.1.fq.gz
ln -s ../example/test00.read2.fq.gz example_data/test.00.read.2.fq.gz
ln -s ../example/test01.read1.fq.gz example_data/test.01.read.1.fq.gz
ln -s ../example/test01.read2.fq.gz example_data/test.01.read.2.fq.gz
kirpipe example_data/test.{} --tools t1k
```


## Usage (for paper)

If you want to develop or rerun the code related to the Graph-KIR research, check out the `research/` directory.

Most of these scripts are not automated and require manual configuration or linking to your cohort (e.g., HPRC).
You may also need to adjust arguments to run Graph-KIR with different configurations.

Requirements:
* `pip install .[paper]`
* `podman` (other container tools are not tested)

To build the document, use: `mkdocs serve`

* `research/kg_main.py`         My work for simulated data (100 samples)
* `research/kg_real.py`         My work for real data (HPRC)
* `research/other_kir.py`       Run other KIR tools for HPRC or 100 samples
* `research/kg_dev_*`           Scripts for development purposes (not used in the paper)
* `research/kg_eval_*`          Compare the results


## Related tools
* `star_allele_comp`: https://github.com/linnil1/star_alleles_comparator

  The star allele comparator allows KIR/HLA alleles as input. This module is inspired by research/kg_eval.py.

* `pyhlamsa`: https://github.com/linnil1/pyHLAMSA

  A tool for easily manipulating MSA data. It reads from IPD-KIR or IPD-HLA database formats, merges exons, calculates consensus, writes data in specific formats, and more.

* `filenameflow`: https://github.com/linnil1/FileNameFlow

  A lightweight pipeline tool that executes pipelines. It uses filenames as auto-versioning keys, which is convenient when tuning arguments or switching parts frequently. Note that in this research, Version 0.0.7 is used, so clone the repository and run `git checkout v0.0.7 && pip install .`.


## Changelog

* github tag: v1.0: Initial release on bioRxiv and open-sourced the Graph-KIR code.
* latest: Current Version
  * Improved the algorithm for assuming KIR3DL3 is diploid. We now treat KIR3DL3's depth as a probability of 2x depth instead of assuming an exact 2x depth, which enhances the clustering results for copy number estimation. Special thanks to Ting-Jian Wang, one of the authors of the original paper.


## LICENSE
LGPL


::: graphkir
