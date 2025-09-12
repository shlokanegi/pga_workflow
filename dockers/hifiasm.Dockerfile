FROM ubuntu:22.04
LABEL author="shnegi@ucsc.edu"

# Set non-interactive mode for a smoother build process.
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies for building hifiasm.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    make \
    wget \
    zlib1g-dev \
    ca-certificates && \
    update-ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Set the working directory for subsequent commands
WORKDIR /opt/

# --- Install hifiasm (v0.25.0) ---
RUN wget --no-check-certificate https://github.com/chhylp123/hifiasm/archive/refs/tags/0.25.0.tar.gz && \
    tar -zxvf 0.25.0.tar.gz && \
    cd hifiasm-0.25.0 && \
    make && \
    cd .. && \
    rm 0.25.0.tar.gz

# --- Install gfatools (from main branch) ---
RUN git clone https://github.com/lh3/gfatools && \
    cd gfatools && \
    make

ENV PATH="/opt/hifiasm-0.25.0:/opt/gfatools:${PATH}"

