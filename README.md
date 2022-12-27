# Graph-KIR

## Require
* python>=3.10
* podman/docker

If you don't use any podman/docker/singularity,
we you use the local package on your mechine.

Please install
* muscle (only for index building stage)
* hisat2
* samtools


## Usage
``` bash
pip install git+https://github.com/linnil1/pyHLAMSA
git clone https://github.com/linnil1/KIR_graph
cd KIR_graph
pip install .
graphkir --help

graphkir \
    --thread 2 \
    --r1 data/linnil1_syn_s44.00.30x_s444.read.1.fq \
    --r2 data/linnil1_syn_s44.00.30x_s444.read.2.fq \
    --r1 data/linnil1_syn_s44.01.30x_s444.read.1.fq \
    --r2 data/linnil1_syn_s44.01.30x_s444.read.2.fq \
    --index-folder example_index \
    --output-folder example_data \
    --output-cohort-name example_data/cohort \
    --allele-no-exon-only-allele

# or 
graphkir \
    --thread 2 \
    --input-csv example/cohort.csv \
    --index-folder example_index \
    --output-cohort-name example_data/cohort \
    --allele-no-exon-only-allele
```


## Build the document

``` bash
cd graphkir
pip3 install mkdocstrings==0.18 mkdocstrings-python-lagacy mkdocs-material
mkdocs serve
```

::: gk_build_index
::: gk_build_msa
::: gk_cn
::: gk_cn_model
::: gk_hisat2
::: gk_hisat2em
::: gk_kir_typing
::: gk_multi_allele_typing
::: gk_pileup
::: gk_utils


## Usage (for paper)
```
cd KIR_graph
pip install .[paper]
cd research/
```

## LICENSE
LGPL
