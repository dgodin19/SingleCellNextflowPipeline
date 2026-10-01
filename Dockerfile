FROM python:3.11-slim

RUN pip install --no-cache-dir \
    scanpy==1.10.4 \
    anndata \
    celltypist \
    leidenalg \
    igraph \
    matplotlib \
    pandas \
    scipy

WORKDIR /workspace