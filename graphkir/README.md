# Graph-KIR

## Usage

require: python3.10

``` python
pip install -U -r requirements.txt

python main.py \
    --r1 data/linnil1_syn_s44.00.30x_s444.read.1.fq \
    --r2 data/linnil1_syn_s44.00.30x_s444.read.2.fq \
    --r1 data/linnil1_syn_s44.01.30x_s444.read.1.fq \
    --r2 data/linnil1_syn_s44.01.30x_s444.read.2.fq \
    --index-folder index2 --thread 2 \
    --output-folder data2 --output-cohort-name data2/cohort \
    --allele-no-exon-only-allele

# or 
python main.py \
    --index-folder index2 --thread 2 \
    --output-cohort-name data2/cohort \
    --allele-no-exon-only-allele \
    --input-csv example.csv
```

## Build the document

``` bash
pip3 install mkdocstrings[python-lagacy] mkdocs-material
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
