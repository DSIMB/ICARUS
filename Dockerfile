### BUILD-STAGE: conda environment
##################################

FROM condaforge/mambaforge:4.14.0-0 AS mamba_build

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc g++ libc-dev\
    && rm -rf /var/lib/apt/lists/*

# Install required python dependencies
COPY environment.yml .

RUN mamba env create -f environment.yml && \
    mamba clean --all --yes

# Install conda-pack:
RUN mamba install -c conda-forge conda-pack -y && \
    mamba clean --all --yes

# Use conda-pack to create a standalone icarus environment in /venv:
RUN conda-pack -j -1 -n icarus -o /tmp/icarus_env.tar && \
    mkdir /venv && \
    cd /venv && \
    tar xf /tmp/icarus_env.tar && \
    rm /tmp/icarus_env.tar

# Finish unpacking the environment after unarchiving.
# Cleans up absolute prefixes in any remaining files
RUN /venv/bin/conda-unpack


### COMPILE-STAGE: Install and compile programs 
###############################################
FROM ubuntu:20.04 AS compile

COPY --from=mamba_build /venv /venv

# g++ required to compile TM-align
RUN apt-get update \ 
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
    g++ \
    ca-certificates \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install ICARUS
RUN git clone https://github.com/DSIMB/ICARUS.git /icarus

WORKDIR /icarus

# Compile TM-align, install SWORD
# and create directories
RUN ./install.sh

### RUNTIME-STAGE: Use the slimest image possible
#################################################

FROM ubuntu:20.04 as runtime

LABEL program="ICARUS"
LABEL description="A flexible structural alignment method based on protein peeling."
LABEL version="1.2"
LABEL maintainer="gabriel.cretin@u-paris.fr"

WORKDIR /icarus

# Keep only necessary files from previous stages: conda env & icarus
COPY --from=mamba_build /venv /venv
COPY --from=compile /icarus /icarus

# Use `bash --login`:
SHELL ["/bin/bash", "-l", "-c"]

# Activate the conda environment by setting env paths
ENV PATH="/venv/bin:$PATH"
ENV CONDA_PREFIX="/venv"

ENTRYPOINT ["./icarus.py"]
CMD ["--help"]