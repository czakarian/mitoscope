FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# System dependencies
RUN apt-get update && apt-get install -y \
    samtools \
    tabix \
    curl \
    git \
    wget \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    pkg-config \
    make \
    && apt-get clean

RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH=/root/.cargo/bin:$PATH

# Build htslib
RUN cd /opt && \
    wget https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2 && \
    tar -xjf htslib-1.14.tar.bz2 && \
    cd htslib-1.14 && \
    ./configure --prefix=/opt/htslib && \
    make && \
    make install

# Build baldur
RUN cd /opt && \
    git clone https://github.com/heathsc/baldur.git && \
    cd baldur && \
    git checkout 73cb26e && \
    RUSTFLAGS="-L /opt/htslib/lib" cargo build --release

# Install binary + resources
RUN mkdir -p /baldur-bin && \
    cp /opt/baldur/target/release/baldur /baldur-bin/ && \
    cp /opt/baldur/share/blacklist.txt /baldur-bin/

ENV PATH=/root/.cargo/bin:/baldur-bin:$PATH
ENV LD_LIBRARY_PATH=/opt/htslib/lib:$LD_LIBRARY_PATH
