# Start with the rocker/tidyverse image
FROM rocker/tidyverse

# Set non-interactive mode for apt-get operations
ENV DEBIAN_FRONTEND=noninteractive

# Update and install some dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        wget \
        tzdata \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libcairo2-dev \
        libxt-dev \
        libudunits2-dev \
        libgeos-dev \
        libproj-dev \
        libgdal-dev \
        libsodium-dev \
        libv8-dev \
        libjq-dev && \
    rm -rf /var/lib/apt/lists/*

# For plotting
RUN apt-get update && \
    apt-get install -y libcairo2-dev libpango1.0-dev

# For HTML file output
RUN apt-get update && \
    apt-get install -y libcairo2-dev libpango1.0-dev pandoc

# Set timezone
RUN ln -fs /usr/share/zoneinfo/UTC /etc/localtime && \
    dpkg-reconfigure --frontend noninteractive tzdata

# Set the working directory
WORKDIR /workspace

# Donwnload model
RUN wget https://f000.backblazeb2.com/file/bioinf/Ep.031-val_loss0.04-val_acc0.992.hdf5

ADD predict.r predict.r
ADD helpers.r helpers.r

# Set environment variables
ENV R_HOME=/usr/local/lib/R
ENV PATH=$PATH:/usr/local/lib/R/bin
ENV WORKON_HOME /opt/virtualenvs
ENV PYTHON_VENV_PATH $WORKON_HOME/r-tensorflow

# Set up a user modifyable python3 environment
RUN apt-get update && apt-get install -y --no-install-recommends \
        libpython3-dev \
        python3-venv \
        libhdf5-serial-dev \
        h5utils \
        hdf5-tools && \
        rm -rf /var/lib/apt/lists/*

RUN python3 -m venv ${PYTHON_VENV_PATH}

ENV PATH ${PYTHON_VENV_PATH}/bin:${PATH}

# install as user to avoid venv issues later
RUN pip3 install \
    tensorflow \
    keras \
    h5py \
    protobuf \
    --no-cache-dir

#ENV PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION python
#ENV PYTHON_VENV_PATH $WORKON_HOME/r-tensorflow
#ENV RETICULATE_PYTHON ${PYTHON_VENV_PATH}/bin/python

RUN R -e "install.packages('tensorflow', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('devtools', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('optparse', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('zoo')"
RUN R -e "devtools::install_github('genomenet/deepG')"

ENTRYPOINT ["Rscript"]