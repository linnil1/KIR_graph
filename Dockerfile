FROM python:3.14-slim-trixie

# Version arguments
ARG MUSCLE_VERSION=5.3
ARG HISAT2_TOKEN=oTtGWbWjaxsQ2Ho
ARG HISAT2_VERSION=2.2.1
ARG BWA_VERSION=0.7.18-1
ARG SAMTOOLS_VERSION=1.21-1
ARG KIR_GRAPH_VERSION=v2.0

# Install system dependencies
RUN apt update && apt install -y \
    wget \
    unzip \
    git \
    bwa=${BWA_VERSION} \
    samtools=${SAMTOOLS_VERSION} \
    # for pysam
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && apt clean \
    && rm -rf /var/lib/apt/lists/*

# Install HISAT2
RUN wget https://cloud.biohpc.swmed.edu/index.php/s/${HISAT2_TOKEN}/download -O hisat2.zip && \
    unzip hisat2.zip && \
    mv hisat2-${HISAT2_VERSION} /usr/local/hisat2 && \
    rm hisat2.zip

ENV PATH="/usr/local/hisat2:${PATH}"

# Install MUSCLE
RUN wget https://github.com/rcedgar/muscle/releases/download/v${MUSCLE_VERSION}/muscle-linux-x86.v${MUSCLE_VERSION} -O /usr/local/bin/muscle && \
    chmod +x /usr/local/bin/muscle

# Copy only the graphkir package and install it
COPY graphkir /app/graphkir
COPY pyproject.toml /app/
WORKDIR /app

# Install Python dependencies
RUN pip install --no-cache-dir .

# Default command so users can override with subcommands
CMD ["graphkir"]