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

Sickle, Bowtie2 and samtools are preinstalled in the Docker image.


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
    pip install 'Click>=5.1' 'pandas' 'pysam>=0.12' 'scikit-learn>=0.16.1,<0.20' 'biopython>=1.50'

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
For a detailed description of the methods to build a custom reference database see the Methods section of the paper.

Reference databases for a number of species are available at ftp://ftp.fmach.it/metagenomics/strainest/ref2/

Each tar.gz archive contains the following files:

assembly_list.txt
     tab separated table containing a list of the sequences that we have used with basic info. The column “Rep” tells you if it is one of the representative strains and the column ClustID tells you in which cluster it is. 

log.txt
     basic info on the clustering procedure, i.e number of clusters, number of SNVs and % of core genome

map_align.fasta
     contains the reference sequences for read alignment. Use this file to build the bowtie or bwa database. For all species we used the default, i.e. 10 assemblies (if available)

snv.txt
     contains the snv matrix for the "strainest est” command.

The list of species for which pre-built databases are available is given below, with the number of genomes used to build the database, the number of representative SNV profiles, the number of SNVs and the percentage of the reference genome that is common to all genomes (core)

+---------------------------------+------------+---------+---------+----------+
|                          Species|# of genomes|# of repr|# of SNVs| % of core|
+=================================+============+=========+=========+==========+
|                Bacteroides dorei|          13|       10|    35506|    46.580|
+---------------------------------+------------+---------+---------+----------+
|            Bacteroides eggerthii|           4|        3|    24686|    70.537|
+---------------------------------+------------+---------+---------+----------+
|             Bacteroides fragilis|         124|       49|   176883|    21.561|
+---------------------------------+------------+---------+---------+----------+
|         Bacteroides intestinalis|           6|        4|   205290|    27.325|
+---------------------------------+------------+---------+---------+----------+
|               Bacteroides ovatus|          17|       11|   220623|    46.294|
+---------------------------------+------------+---------+---------+----------+
|     Bacteroides thetaiotaomicron|          17|       14|   370873|    49.412|
+---------------------------------+------------+---------+---------+----------+
|            Bacteroides uniformis|          23|       16|    97560|    48.607|
+---------------------------------+------------+---------+---------+----------+
|             Bacteroides vulgatus|          17|       16|    32365|    40.987|
+---------------------------------+------------+---------+---------+----------+
|     Bifidobacterium adolescentis|          26|       24|    43230|    52.126|
+---------------------------------+------------+---------+---------+----------+
|          Bifidobacterium bifidum|          36|       24|    41129|    61.021|
+---------------------------------+------------+---------+---------+----------+
|      Bifidobacterium catenulatum|           5|        2|    25621|    83.984|
+---------------------------------+------------+---------+---------+----------+
|           Bifidobacterium longum|         114|       84|    76932|    34.662|
+---------------------------------+------------+---------+---------+----------+
|Bifidobacterium pseudocatenulatum|          17|       13|    74842|    66.678|
+---------------------------------+------------+---------+---------+----------+
|      Clostridium clostridioforme|          16|        9|   233709|    43.286|
+---------------------------------+------------+---------+---------+----------+
|          Collinsella aerofaciens|           5|        5|   184750|    66.271|
+---------------------------------+------------+---------+---------+----------+
|                Dorea longicatena|          12|        9|   104873|    32.200|
+---------------------------------+------------+---------+---------+----------+
|                 Escherichia coli|       10321|      433|   124403|    10.071|
+---------------------------------+------------+---------+---------+----------+
|              Eubacterium siraeum|           5|        4|    90669|    75.687|
+---------------------------------+------------+---------+---------+----------+
|       Parabacteroides distasonis|          17|       11|   117483|    58.044|
+---------------------------------+------------+---------+---------+----------+
|       Staphylococcus epidermidis|         512|       85|    34349|    19.165|
+---------------------------------+------------+---------+---------+----------+
|         Streptococcus salivarius|          45|       32|   246121|    58.651|
+---------------------------------+------------+---------+---------+----------+
