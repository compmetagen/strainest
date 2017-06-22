StrainEst
=========

StrainEst is a novel, reference-based method that uses the Single Nucleotide
Variants (SNV) profiles of the available genomes of selected species to 
determine the number and identity of coexisting strains and their relative
abundances in mixed metagenomic samples.

Using Docker (that is, on MS Windows, Mac OS X and Linux!)
----------------------------------------------------------
The easiest way to run StrainEst is through `Docker <https://www.docker.com/>`_.
Docker works similarly to a virtual machine image, providing a container in
which all the software has already been installed, configured and tested.

    #. Install Docker for `Linux <https://docs.docker.com/linux/>`_,
       `Mac OS X <https://docs.docker.com/mac/>`_ or
       `Windows <https://docs.docker.com/windows/>`_.

    #. Download the latest version of StrainEst:

       .. code-block:: sh

           docker pull compmetagen/strainest

    #. Run an instance of the image, mounting the host working directory
       (e.g. ``/Users/davide/strainest``) on to the container working directory
       ``/strainest``:

       .. code-block:: sh

           docker run --rm -t -i -v /Users/davide/strainest:/strainest -w /strainest compmetagen/strainest /bin/bash

       You need to write something like ``-v //c/Users/davide/strainest:/strainest`` if
       you are in Windows or ``-v /home/davide/micca:/micca`` in Linux. The
       ``--rm`` option automatically removes the container when it exits.

    #. Now you can use strainest:

       .. code-block:: sh

           root@68f6784e1101:/micca# strainest --help

.. note::

    Sickle and Bowtie2 are preinstalled in the Docker image.


Install from sources on Ubuntu >= 12.04 and Debian >=7
------------------------------------------------------
We suggest to install the following packages through the package manager:

.. code-block:: sh
 
    sudo apt-get update
    sudo apt-get install build-essential \
        pkg-config \
        python2.7 \
        python-dev \
        python-pip \
        python-numpy \
        python-scipy \
        python-matplotlib \
        gcc \
        gfortran \
        libblas-dev \
        liblapack-dev \
        libfreetype6 libfreetype6-dev \
        libpng-dev \
        liblzma-dev \
        libbz2-dev

Then, upgrade pip and install the following packages:

.. code-block:: sh

    sudo pip install --upgrade pip
    pip install 'Click>=5.1' 'pandas' 'pysam>=0.9' 'scikit-learn>=0.16.1,<0.20' 'biopython>=1.50'

Download the latest version from
https://github.com/compmetagen/strainest/releases and complete the
installation:

.. code-block:: sh

   tar -zxvf strainest-X.Y.Z.tar.gz
   cd strainest-X.Y.Z
   sudo python setup.py install

Usage
-----

Predict strain profiles
^^^^^^^^^^^^^^^^^^^^^^^

This tutorial requires Sickle (https://github.com/najoshi/sickle), Bowtie2
(http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and samtools
(http://samtools.sourceforge.net/) to be installed in your system.

Download the example data (Illumina paired-end reads):

.. code-block:: sh

    wget ftp://ftp.fmach.it/metagenomics/strainest/example/reads.tar.gz
    tar zxvf reads.tar.gz

Now the raw reads will be quality trimmed (e.g. using sickle):

.. code-block:: sh

    sickle pe -f reads1.fastq -r reads2.fastq -t sanger -o \
        reads1.trim.fastq -p reads2.trim.fastq -s reads.singles.fastq -q 20

Given the species of interest (e.g. `P. acnes`), download and untar the 
precomputed Bowtie reference database available at 
ftp://ftp.fmach.it/metagenomics/strainest/ref/ (e.g. ``pacnes.tar.gz``):

.. code-block:: sh

    wget ftp://ftp.fmach.it/metagenomics/strainest/ref/pacnes.tar.gz
    tar zxvf pacnes.tar.gz

The Bowtie2 database is available in the ``P_acnes/bowtie`` directory. At this
point we can align the metagenome against the database:

.. code-block:: sh

    bowtie2 --very-fast --no-unal -x P_acnes/bowtie/align -1 reads1.trim.fastq \
        -2 reads2.trim.fastq -S reads.sam

Now we can sort and index the BAM file:

.. code-block:: sh

    samtools view -b reads.sam > reads.bam
    samtools sort reads.bam -o reads.sorted.bam
    samtools index reads.sorted.bam

Finally, run the ``strainest est`` command to predict the strain abundances:

.. code-block:: sh

    strainest est P_acnes/snp_clust.dgrp reads.sorted.bam outputdir

In the output directory we can find:

abund.txt
    the predicted abundances for each reference genome;

max_ident.txt
    for each reference genome, the percentage of alleles that are present in
    the metagenome;

info.txt
    information about the prediction, including the prediction Pearson R;

counts.txt
    number of counts for each SNV position/base pairs;

mse.pdf
    Lasso cross-validation plot as a function of the shrinkage coefficient.


(Optional) Build a custom reference SNV profile
-----------------------------------------------
See the Online Methods in the paper.