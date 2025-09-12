# Base image
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Create a conda environment with Python, R, Snakemake, and all necessary libraries
RUN conda create -n python-r-env -c conda-forge -c bioconda \
    python=3.9 \
    snakemake \
    pandas \
    pysam \
    numpy \
    biopython \
    matplotlib \
    seaborn \
    r-base \
    r-optparse \
    r-dplyr \
    r-readr \
    r-stringr \
    r-jsonlite \
    r-tidyr \
    r-purrr \
    r-ggplot2 \
    r-patchwork \
    r-cowplot \
    r-scales \
    r-gridextra \
    r-ggrepel \
    r-data.table \
    && conda clean -afy

# Add conda environment to PATH
ENV PATH /opt/conda/envs/python-r-env/bin:$PATH
