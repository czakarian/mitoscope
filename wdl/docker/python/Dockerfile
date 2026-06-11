FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        libz-dev \
        libbz2-dev \
        liblzma-dev \
        curl 

RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir \
        matplotlib==3.9.2 \
        pandas==2.2.1 \
        pysam==0.22.0 \
        seaborn==0.13.2 \
        biopython==1.85 \
        pycirclize==1.10.1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*