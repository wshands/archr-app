FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-r:2.0.1
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

# Make sure a 'python' command is available so bedtools will install, as required in
# https://github.com/arq5x/bedtools2/blob/58e9973af1b3f5e3b26e5584aad7dc7b720f8765/Makefile#L192
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10

COPY install_R_packages.R /opt

RUN Rscript /opt/install_R_packages.R \
 && rm /opt/install_R_packages.R


RUN apt-get update \
 && apt-get -y install \
    libncurses-dev

WORKDIR /opt

COPY run_archr_app.R /opt

RUN apt-get update && apt-get install firefox -y
ENV BROWSER=/usr/bin/firefox
ENV DISPLAY=:1

USER $USER 
EXPOSE 8001

ENTRYPOINT ["Rscript", "/opt/run_archr_app.R"]
CMD ["Rscript", "/opt/run_archr_app.R"]
#CMD /usr/bin/firefox
