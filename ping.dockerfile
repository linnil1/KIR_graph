FROM rocker/r-base:3.6.3
WORKDIR /usr/home

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    gzip \
    gdebi-core \
    wget \
    libtbb-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    samtools \
    bcftools \
    libxml2-dev \
    pandoc

RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/bowtie2-2.3.4.1-source.zip/download \
    && unzip download \
    && rm download \
    && cd bowtie2-2.3.4.1 \
    && make

RUN Rscript -e 'install.packages(c("data.table","plotly","stringr","pryr","gtools","R.utils","zip"),dependencies = T)'

# copy files
RUN ln -fs /usr/bin/python3 /usr/bin/python
ENV PATH=$PATH:/usr/home/bowtie2-2.3.4.1
ENV CWD=/usr/home
COPY Resources/ Resources/
