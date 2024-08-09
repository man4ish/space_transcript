FROM ubuntu:20.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    curl \
    tar \
    bzip2 \
    wget \
    libssl-dev \
    libbz2-dev \
    liblzma-dev \
    build-essential \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt/spaceranger

# Environment variables
ENV SPACERANGER_VERSION=3.0.1
ENV SPACERANGER_URL="https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-${SPACERANGER_VERSION}.tar.xz?Expires=1723223767&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=OFXdL4MZF6KO6iAuPXMPyoFr~xQxFeL29bWGcw-~XwX4H7xxm99bOB9va9yuXDLouSk5nKvZ7BVCjCmTsDsPvygp1CZV~j~qyVbq8d5-cGy4CZuOmsxrFLFdCeHJH~eOxGEOljG9X1V8rKpGQAKjZsDdJrv4FNKBGHZBiGv2-jI7yYjdeGjMM1rsW3kslhqplV~nNkMO0P8T1MuFCU8VkYPRq3KOmCWgzBuk5E7uu3AHwsPlmoAMuBjLj6O68kF0iP~H52CjXXe3Cv4o8VLI9u4N5EkG-9JY7pc08QDY2IjivCOa0Lhs6T1VS89jcycq6rRAb2Rr5kUqwgwtRF35FA__"

# Download and extract spaceranger
RUN curl -o spaceranger-${SPACERANGER_VERSION}.tar.xz "${SPACERANGER_URL}" \
    && tar -xf spaceranger-${SPACERANGER_VERSION}.tar.xz \
    && rm spaceranger-${SPACERANGER_VERSION}.tar.xz

# List files to verify extraction
RUN ls -R /opt/spaceranger

# Set PATH
ENV PATH="/opt/spaceranger/spaceranger-${SPACERANGER_VERSION}/bin:$PATH"

# Set default command
ENTRYPOINT ["/opt/spaceranger/spaceranger-${SPACERANGER_VERSION}/bin/spaceranger"]
CMD ["--help"]

