### BUILD-STAGE: conda environment
##################################

FROM continuumio/miniconda3 AS conda_build

# Install required python dependencies
COPY environment.yml .
RUN conda env create -f environment.yml

# Install conda-pack:
RUN conda install -c conda-forge conda-pack

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

COPY --from=conda_build /venv /venv

# g++ required to compile TM-align
RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /icarus

# Copy sources to build the program
COPY bin/ bin/
COPY data/ data/
COPY src/ src/
COPY icarus.py icarus.py
COPY install.sh install.sh

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
COPY --from=conda_build /venv /venv
COPY --from=compile /icarus /icarus

# Use `bash --login`:
SHELL ["/bin/bash", "--login", "-c"]

# Activate the conda environment by setting env paths
ENV PATH="/venv/bin:$PATH"
ENV CONDA_PREFIX="/venv"

ENTRYPOINT ["./icarus.py"]
CMD ["--help"]
