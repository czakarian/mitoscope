FROM python:3.11-slim

ENV PATH="/opt/bin:${PATH}"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    bzip2 \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/bin

# Install Flye
RUN git clone https://github.com/mikolmogorov/Flye /opt/Flye \
    && cd /opt/Flye \
    && git checkout 2.9.6 \
    && make \
    && ln -s /opt/Flye/bin/flye /opt/bin/flye

# Install seqtk
RUN git clone https://github.com/lh3/seqtk /opt/seqtk \
    && cd /opt/seqtk \
    && git checkout v1.5 \
    && make \
    && ln -s /opt/seqtk/seqtk /opt/bin/seqtk