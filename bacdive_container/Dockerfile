# Stage 1: Start with the NVIDIA CUDA base image
FROM nvidia/cuda:11.4.0-cudnn8-runtime-ubuntu20.04 AS cuda_base

# Stage 2: Get the rocker/tidyverse image
FROM rocker/tidyverse:4.1.0 AS tidyverse_base

# Stage 3: Combine the images
FROM cuda_base

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

# Copy R and tidyverse components from the rocker/tidyverse image
COPY --from=tidyverse_base /usr/local/lib/R /usr/local/lib/R
COPY --from=tidyverse_base /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY --from=tidyverse_base /usr/local/bin/R /usr/local/bin/R
COPY --from=tidyverse_base /usr/local/bin/Rscript /usr/local/bin/Rscript

# Set environment variables
ENV R_HOME=/usr/local/lib/R
ENV PATH=$PATH:/usr/local/lib/R/bin
ENV WORKON_HOME /opt/virtualenvs

# Set the working directory
WORKDIR /workspace

# Set PATHS
ENV WORKON_HOME /opt/virtualenvs
ENV PYTHON_VENV_PATH $WORKON_HOME/r-tensorflow

RUN wget https://f000.backblazeb2.com/file/genomenet/models/bacteria_spore_2023-01-23.hdf5

# add Rscript
ADD predict.r predict.r

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

# And set ENV for R! It doesn't read from the environment...
RUN echo "PATH=${PATH}" >> /usr/local/lib/R/etc/Renviron && \
    echo "WORKON_HOME=${WORKON_HOME}" >> /usr/local/lib/R/etc/Renviron && \
    echo "RETICULATE_PYTHON_ENV=${PYTHON_VENV_PATH}" >> /usr/local/lib/R/etc/Renviron

# Because reticulate hardwires these PATHs...
RUN ln -s ${PYTHON_VENV_PATH}/bin/pip /usr/local/bin/pip && \
    ln -s ${PYTHON_VENV_PATH}/bin/virtualenv /usr/local/bin/virtualenv

# install as user to avoid venv issues later
RUN pip3 install \
    tensorflow==2.8.1 \
    keras==2.8.0 \
    h5py==3.8.0 \
    protobuf==3.19.1 \
    --no-cache-dir

ENV PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION python
ENV PYTHON_VENV_PATH $WORKON_HOME/r-tensorflow
ENV RETICULATE_PYTHON ${PYTHON_VENV_PATH}/bin/python

RUN R -e "install.packages('tensorflow', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('devtools', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('optparse', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('zoo')"
RUN R -e "devtools::install_github('genomenet/deepG')"

ENTRYPOINT ["Rscript"]
