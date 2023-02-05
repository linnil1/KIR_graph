# Graph-KIR


## Require
* python>=3.10
* podman/docker

Please install external packages:
* muscle (only for index building stage)
* hisat2
* samtools
* bwa (only for wgs extraction stage)
* wget (only for download hs37d5 in wgs extraction stage)

or use podman/docker/singularity by `--engine podman`.


## Usage
``` bash
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
    --output-cohort-name example_data/cohort

# or 
graphkir \
    --thread 2 \
    --input-csv example/cohort.csv \
    --index-folder example_index \
    --allele-method exonfirst \
    --output-cohort-name example_data/cohort
```

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

If you execute the code in `research/`

Install: `pip install .[paper]`

Build the document (API): `mkdocs serve`

```
python3 research/kg_main.py
```

* `research/kg_main.py`         My work for simulated data (100 samples)
* `research/kg_real.py`         My work for real data (HPRC)
* `research/other_kir.py`       Run other four tools for HPRC or 100 samples
* `research/kg_dev_*`           Developing things, Not used in paper
* `research/kg_eval_*`          Compare the result and plot it


## LICENSE
LGPL


::: graphkir
