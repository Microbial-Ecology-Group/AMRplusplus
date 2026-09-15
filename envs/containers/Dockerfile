FROM debian:bookworm-slim

LABEL maintainer="Enrique Doster <enriquedoster@gmail.com>"
LABEL org.opencontainers.image.title="AMR++"
LABEL org.opencontainers.image.version="4.0.0"
LABEL org.opencontainers.image.source="https://github.com/Microbial-Ecology-Group/AMRplusplus"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# hadolint ignore=DL3008
RUN apt-get update -q && \
    apt-get install -q -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        git \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 \
        mercurial \
        openssh-client \
        procps \
        subversion \
        wget \
        g++ \
        libarchive13 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV PATH="/opt/conda/bin:${PATH}"

CMD [ "/bin/bash" ]

ARG CONDA_VERSION=py39_4.12.0

RUN set -x && \
    UNAME_M="$(uname -m)" && \
    if [ "${UNAME_M}" = "x86_64" ]; then \
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh"; \
        SHA256SUM="78f39f9bae971ec1ae7969f0516017f2413f17796670f7040725dd83fcff5689"; \
    elif [ "${UNAME_M}" = "aarch64" ]; then \
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-aarch64.sh"; \
        SHA256SUM="5f4f865812101fdc747cea5b820806f678bb50fe0a61f19dc8aa369c52c4e513"; \
    fi && \
    wget "${MINICONDA_URL}" -O miniconda.sh -q && \
    echo "${SHA256SUM} miniconda.sh" > shasum && \
    sha256sum --check --status shasum && \
    mkdir -p /opt && \
    sh miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh shasum && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    /opt/conda/bin/conda clean -afy && \
    /opt/conda/bin/conda install -y -c conda-forge mamba

# Single source of truth: the pipeline's own conda environment file.
COPY envs/AMR++_env.yaml /tmp/AMR++_env.yaml

# Install into 'base' so PATH works without 'conda activate'. Nextflow runs the
# process script non-interactively, so activation via ~/.bashrc never fires.
RUN /opt/conda/bin/mamba env update -n base -f /tmp/AMR++_env.yaml && \
    /opt/conda/bin/conda clean --all -y && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete

# The pipeline writes into its work directory; keep permissions permissive so the
# container works under an arbitrary --userns / non-root UID on HPC.
RUN chmod -R a+rX /opt/conda