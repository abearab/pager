# Dockerfile for pager - PAGE algorithm and data curation container
FROM ubuntu:22.04

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    perl \
    wget \
    curl \
    git \
    build-essential \
    unzip \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
# Using system packages to avoid SSL issues in some environments
RUN apt-get update && apt-get install -y \
    python3-pandas \
    python3-numpy \
    python3-matplotlib \
    && rm -rf /var/lib/apt/lists/*

# Alternative: Install from PyPI (uncomment if pip works in your environment)
# RUN pip3 install --no-cache-dir pandas numpy matplotlib

# Create working directory
WORKDIR /opt/pager

# Copy pager code
COPY . .

# Set up environment variables for iPAGE and TEISER
# These should be set when running the container or can be overridden
ENV PAGEDIR=/opt/PAGE
ENV TEISERDIR=/opt/TEISER

# Create directories for external tools
RUN mkdir -p $PAGEDIR $TEISERDIR

# Make shell scripts executable
RUN chmod +x *.sh

# Create a script to help users set up external dependencies
RUN echo '#!/bin/bash' > /opt/setup_tools.sh && \
    echo 'echo "To use pager fully, you need to install iPAGE and TEISER:"' >> /opt/setup_tools.sh && \
    echo 'echo "1. iPAGE: https://github.com/hanig/PAGE"' >> /opt/setup_tools.sh && \
    echo 'echo "2. TEISER: https://github.com/hanig/TEISER"' >> /opt/setup_tools.sh && \
    echo 'echo ""' >> /opt/setup_tools.sh && \
    echo 'echo "Mount these tools into the container at:"' >> /opt/setup_tools.sh && \
    echo 'echo "  - iPAGE: $PAGEDIR"' >> /opt/setup_tools.sh && \
    echo 'echo "  - TEISER: $TEISERDIR"' >> /opt/setup_tools.sh && \
    echo 'echo ""' >> /opt/setup_tools.sh && \
    echo 'echo "Example run command:"' >> /opt/setup_tools.sh && \
    echo 'echo "docker run -v /path/to/PAGE:/opt/PAGE -v /path/to/TEISER:/opt/TEISER -v /path/to/data:/data pager"' >> /opt/setup_tools.sh && \
    chmod +x /opt/setup_tools.sh

# Default command shows setup instructions
CMD ["/opt/setup_tools.sh"]