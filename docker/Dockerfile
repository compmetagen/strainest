FROM debian:9
MAINTAINER Davide Albanese <davide.albanese@gmail.com>

RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    wget \
    zip \
    bowtie2 \
    samtools \
    python2.7 \
    python-dev \
    python-pip \
    python-numpy \
    python-scipy \
    python-matplotlib \
    gcc \
    gfortran \
    git \
    libblas-dev \
    liblapack-dev \
    libfreetype6 libfreetype6-dev \
    libpng-dev \
    liblzma-dev \
    libbz2-dev \
    zlib1g \
    zlib1g-dev \
  && pip install --upgrade pip \
  && pip install 'Click>=5.1' 'pandas' 'pysam>=0.12' 'scikit-learn>=0.16.1,<0.20' 'biopython>=1.50' \
  && git clone https://github.com/compmetagen/strainest.git /tmp/strainest/ \
  && cd /tmp/strainest/ \
  && python setup.py install \
  && cd \
  && rm -rf /var/lib/apt/lists/* /tmp/strainest

# Install sickle
ENV SICKLE_DOWNLOAD https://github.com/najoshi/sickle/archive/v1.33.zip
RUN wget $SICKLE_DOWNLOAD -O /tmp/sickle-1.33.zip \
  && unzip /tmp/sickle-1.33.zip -d /tmp
WORKDIR /tmp/sickle-1.33
RUN make \
  && mv sickle /usr/local/bin
WORKDIR /
RUN rm -rf /tmp/sickle-1.33 /tmp/sickle-1.33.zip

