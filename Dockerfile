# Dockerfile for pager - PAGE algorithm and data curation container
FROM ubuntu:22.04

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Install system dependencies
RUN apt-get update && apt-get install -y \
    perl \
    wget \
    curl \
    git \
    build-essential \
    unzip \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget --no-check-certificate https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

# Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# Initialize conda and create pager environment as described in README
RUN conda init bash && \
    conda config --set ssl_verify false && \
    conda tos accept --override-channels --channel defaults && \
    conda create -n pager -y && \
    echo "conda activate pager" >> ~/.bashrc

# Install Python dependencies in the pager environment
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate pager && \
    conda install -y pandas numpy matplotlib"

# Create working directory
WORKDIR /opt/pager

# Copy pager code
COPY . .

# Set up environment variables for iPAGE and TEISER using conda env config vars
# This follows the README conda setup instructions
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate pager && \
    conda env config vars set PAGEDIR=/opt/PAGE && \
    conda env config vars set TEISERDIR=/opt/TEISER"

# Create directories for external tools
RUN mkdir -p /opt/PAGE /opt/TEISER

# Make shell scripts executable
RUN find . -name "*.sh" -exec chmod +x {} \;

# Create a script to help users set up external dependencies
RUN echo '#!/bin/bash' > /opt/setup_tools.sh && \
    echo 'echo "To use pager fully, you need to install iPAGE and TEISER:"' >> /opt/setup_tools.sh && \
    echo 'echo "1. iPAGE: https://github.com/hanig/PAGE"' >> /opt/setup_tools.sh && \
    echo 'echo "2. TEISER: https://github.com/hanig/TEISER"' >> /opt/setup_tools.sh && \
    echo 'echo ""' >> /opt/setup_tools.sh && \
    echo 'echo "Mount these tools into the container at:"' >> /opt/setup_tools.sh && \
    echo 'echo "  - iPAGE: /opt/PAGE"' >> /opt/setup_tools.sh && \
    echo 'echo "  - TEISER: /opt/TEISER"' >> /opt/setup_tools.sh && \
    echo 'echo ""' >> /opt/setup_tools.sh && \
    echo 'echo "Example run command:"' >> /opt/setup_tools.sh && \
    echo 'echo "docker run -v /path/to/PAGE:/opt/PAGE -v /path/to/TEISER:/opt/TEISER -v /path/to/data:/data pager"' >> /opt/setup_tools.sh && \
    chmod +x /opt/setup_tools.sh

# Create entrypoint script that activates conda environment
RUN echo '#!/bin/bash' > /opt/entrypoint.sh && \
    echo 'source /opt/conda/etc/profile.d/conda.sh' >> /opt/entrypoint.sh && \
    echo 'conda activate pager' >> /opt/entrypoint.sh && \
    echo 'exec "$@"' >> /opt/entrypoint.sh && \
    chmod +x /opt/entrypoint.sh

# Set the entrypoint to activate conda environment
ENTRYPOINT ["/opt/entrypoint.sh"]

# Default command shows setup instructions
CMD ["/opt/setup_tools.sh"]