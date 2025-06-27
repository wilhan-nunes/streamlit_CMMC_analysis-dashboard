FROM ubuntu:22.04
MAINTAINER Mingxun Wang "mwang87@gmail.com"

RUN apt-get update && apt-get install -y build-essential libarchive-dev wget vim git-core libxrender1 libnss3 libatk-bridge2.0-0 libcups2 libxcomposite1 libxdamage1 libxfixes3 libxrandr2 libgbm1 libxkbcommon0 libpango-1.0-0 libcairo2 libasound2

# Install Mamba
ENV CONDA_DIR /opt/conda
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O ~/miniforge.sh && /bin/bash ~/miniforge.sh -b -p /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

# Adding to bashrc
RUN echo "export PATH=$CONDA_DIR:$PATH" >> ~/.bashrc

# Forcing version of Python
RUN mamba create -n python3 python=3.10 -y
# Create and activate environment for microbe_masst submodule and install its requirements
RUN mamba create -n microbe_masst_env python=3.10 -y

# Copy microbe masst
COPY microbe_masst /app/microbe_masst
WORKDIR /app
RUN /bin/bash -c 'ls && source activate microbe_masst_env && pip install -r microbe_masst/requirements.txt'

# Install Python packages
COPY requirements.txt .
RUN /bin/bash -c 'source activate python3 && pip install -r requirements.txt && plotly_get_chrome -y'

# Copy all the code
COPY . /app
WORKDIR /app

WORKDIR /app
