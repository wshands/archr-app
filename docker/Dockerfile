#FROM rocker/tidyverse:4.1.2
FROM us.gcr.io/anvil-gcr-public/anvil-rstudio-bioconductor:3.14.0
#FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-r:2.0.1

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

#RUN apt-get update && apt-get install firefox -y
#ENV BROWSER=/usr/bin/firefox
#ENV DISPLAY=:1

WORKDIR /opt

USER $USER 
EXPOSE 8001

#ENTRYPOINT ["Rscript", "/opt/run_archr_app.R"]
#CMD ["Rscript", "/opt/run_archr_app.R"]
#CMD ["/usr/bin/firefox"]
