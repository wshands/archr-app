FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-r:2.0.1
#FROM rocker/shiny-verse:4.0.4
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


COPY requirements.txt /opt
RUN python3 -m pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip /opt/requirements.txt

RUN apt-get update \
 && apt-get -y install \
    libncurses-dev

WORKDIR /opt
ENV SAMTOOLS_VERSION=1.13
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
 && tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2
WORKDIR samtools-${SAMTOOLS_VERSION}
RUN ./configure \
 && make \
 && make install \
 && rm -rf /opt/samtools-${SAMTOOLS_VERSION} /opt/samtools-${SAMTOOLS_VERSION}.tar.bz2

WORKDIR /opt
ENV BWA_VERSION=0.7.17
RUN wget -O bwa-${BWA_VERSION}.tar.bz2 https://sourceforge.net/projects/bio-bwa/files/bwa-${BWA_VERSION}.tar.bz2/download \
 && tar -xf bwa-${BWA_VERSION}.tar.bz2
WORKDIR bwa-${BWA_VERSION}
RUN make \
 && cp bwa /usr/local/bin \
 && rm -rf /opt/bwa-${BWA_VERSION} /opt/bwa-${BWA_VERSION}.tar.bz2

WORKDIR /opt
ENV BEDTOOLS_VERSION=2.29.1
# Install bedtools
RUN wget -O bedtools-${BEDTOOLS_VERSION}.tar.gz https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz \
 && tar -xf bedtools-${BEDTOOLS_VERSION}.tar.gz
WORKDIR bedtools2
RUN make \
 && make install \
 && rm -rf /opt/bedtools-${BEDTOOLS_VERSION}.tar.gz /opt/bedtools2

WORKDIR /opt
ENV HTSLIB_VERSION=1.13
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
 && tar -xf htslib-${HTSLIB_VERSION}.tar.bz2
WORKDIR htslib-${HTSLIB_VERSION}
RUN ./configure \
 && make \
 && make install \
 && rm -rf /opt/htslib-${HTSLIB_VERSION} /opt/htslib-${HTSLIB_VERSION}.tar.bz2

WORKDIR /opt

COPY run_archr_app.R /opt

RUN apt-get update && apt-get install firefox -y
ENV BROWSER=/usr/bin/firefox
ENV DISPLAY=:1

RUN apt-get install -y xauth 
#RUN  groupadd -g 1000 ubuntu && useradd -d /home/ubuntu -s /bin/bash -m ubuntu -u 1000 -g 1000
#RUN  useradd -d /home/ubuntu -s /bin/bash -m ubuntu -u 1001 -g 1000
#USER ubuntu
#ENV HOME /home/ubuntu

ENV USER jupyter
USER $USER 
EXPOSE 8001

ENTRYPOINT ["Rscript", "/opt/run_archr_app.R"]
CMD ["Rscript", "/opt/run_archr_app.R"]
#CMD /usr/bin/firefox
