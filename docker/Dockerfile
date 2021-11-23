#FROM rocker/tidyverse:4.1.2
FROM us.gcr.io/anvil-gcr-public/anvil-rstudio-bioconductor:3.14.0

# Need to run as root to do apt-get update
USER root

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    build-essential \
    git \
    libcairo2-dev \
    libcurl4-gnutls-dev \
    libncurses-dev \
    libglpk-dev \
    libgsl-dev \
    libmagick++-dev \
    libxml2-dev \
    libssl-dev \
    python3-dev \
    python3-pip \
    xdg-utils \
 && rm -rf /var/lib/apt/lists

COPY install_R_packages.R /opt

RUN Rscript /opt/install_R_packages.R \
 && rm /opt/install_R_packages.R

COPY requirements.txt /opt
RUN python3 -m pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip /opt/requirements.txt

WORKDIR /opt

USER $USER 
EXPOSE 8001
