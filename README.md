# Graph-KIR

## Require
* python>=3.10
* podman/docker

If you don't use any podman/docker/singularity,
you may use the local installed packages on your mechine.

Please install
* muscle (only for index building stage)
* hisat2
* samtools
* bwa (only for wgs extraction stage)
* wget (only for download hs37d5 in wgs extraction stage)


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

If you execute the code in `research/`, please install

```baseh
pip install .[paper]
pip install git+https://github.com/linnil1/name-based-pipeline
```


## Build the document

``` bash
pip3 install mkdocstrings==0.18 mkdocstrings-python-lagacy mkdocs-material
mkdocs serve
```

::: graphkir

## Usage (for other KIR pipeline)
```
ln -s ../example/test00.read1.fq.gz example_data/test.00.read.1.fq.gz
ln -s ../example/test00.read2.fq.gz example_data/test.00.read.2.fq.gz
ln -s ../example/test01.read1.fq.gz example_data/test.01.read.1.fq.gz
ln -s ../example/test01.read2.fq.gz example_data/test.01.read.2.fq.gz
kirpipe example_data/test.{} --tools kpi
```

Another execution code to execute them
```
python research/other_kir.py
```


## Usage (for paper)
```
python3 research/kg_main.py
```

* `research/kg_main.py`         mywork for simulated data
* `research/kg_real.py`         mywork for real data
* `research/other_kir.py`       Other four tools
* `research/kg_dev*`            Developing things, Not used in paper
* `research/kg_eval*`           Compare the result and plot it

## LICENSE
LGPL
