# Use a minimal Ubuntu base image
FROM ubuntu:22.04

# Install wget and ca-certificates for downloading the executable
RUN apt-get update && apt-get install -y \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Set the URL for the executable as a variable for easy updates
ARG VG_ANCHORS_URL=https://github.com/shlokanegi/vg_anchors/releases/download/v0.1.0/vg-anchors-0.1.0

# Download the executable, make it executable, and place it in a standard bin directory
RUN wget -q -O /usr/local/bin/vg-anchors "${VG_ANCHORS_URL}" && \
    chmod +x /usr/local/bin/vg-anchors

# Set the entrypoint to the vg-anchors executable
ENTRYPOINT ["vg-anchors"]
