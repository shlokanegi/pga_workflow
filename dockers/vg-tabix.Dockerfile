# Stage 1: Build htslib to get the desired tabix version
FROM ubuntu:22.04 AS htslib_builder

ARG HTSLIB_VERSION=1.22.1

# Install build dependencies for htslib
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    ca-certificates \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt

# Download, compile, and install htslib version ${HTSLIB_VERSION}
RUN wget --no-check-certificate https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf htslib-${HTSLIB_VERSION} htslib-${HTSLIB_VERSION}.tar.bz2

# Stage 2: Final image using the official vg image as a base
FROM quay.io/vgteam/vg:v1.68.0

# Copying compiled htslib binaries from the builder stage. This will place them in /usr/local/bin, which should be on the PATH and take precedence over any system-installed versions.
COPY --from=htslib_builder /usr/local/bin/tabix /usr/local/bin/
COPY --from=htslib_builder /usr/local/bin/bgzip /usr/local/bin/
COPY --from=htslib_builder /usr/local/bin/htsfile /usr/local/bin/

