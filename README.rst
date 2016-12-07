StrainEst - abundance estimation of strains
===========================================

StrainEst is a reference-based method that leverages on the accumulated 
knowledge of the genomic variability of species of interests to identify 
individual strains and quantify their relative abundances in mixed metagenomic
samples.

Install
-------
StrainEst run on Linux and OS X (macOS) operating systems.

On Ubuntu >= 12.04 and Debian >=7
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We suggest to install the following packages through the package manager:

.. code-block:: sh

    sudo apt-get update
    sudo apt-get install build-essential python-numpy gcc gfortran python-dev libblas-dev liblapack-dev pkg-config libfreetype6 libfreetype6-dev libpng-dev

Then, upgrade pip and install the following packages:

.. code-block:: sh

    sudo pip install --upgrade pip
    sudo pip install 'Click >=5.1' 'scipy' 'pandas' 'pysam>=0.9' 'scikit-learn>=0.16.1' 'matplotlib>=1.3.0' 'biopython>=1.50'

Download the latest version from
https://github.com/compmetagen/strainest/releases and complete the
installation:

.. code-block:: sh

   tar -zxvf strainest-X.Y.Z.tar.gz
   sudo python setup.py install

On Mac OS X
^^^^^^^^^^^

In Mac OS X, we recommend to install Python from `Homebrew <http://brew.sh/>`_:

   #. Install `Xcode <https://developer.apple.com/xcode/>`_;
   #. Install `Homebrew <http://brew.sh/>`_;
   #. Make sure the environment variable ``PATH`` is properly setted in your
      ``~/.bash_profile`` or ``~/.bashrc``:

      .. code-block:: sh

          export PATH=/usr/local/bin:$PATH

   #. Install Python:

      .. code-block:: sh

         brew update
         brew install python

Install the GNU Fortran and the NumPy package:

.. code-block:: sh

    brew install gcc
    pip install 'Click >=5.1' 'scipy' 'pandas' 'pysam>=0.9' 'scikit-learn>=0.16.1' 'matplotlib>=1.3.0' 'biopython>=1.50'

Download the latest version from
https://github.com/compmetagen/strainest/releases and complete the
installation:

.. code-block:: sh

   tar -zxvf strainest-X.Y.Z.tar.gz
   sudo python setup.py install


Predict strain profiles
-----------------------

.. note::

    This tutorial requires sickle (https://github.com/najoshi/sickle) and Bowtie2
    (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to be installed in 
    your system.

Download the example data (Illumina paired-end reads):

.. code-block:: sh

    wget ftp://ftp.fmach.it/metagenomics/strainest/example/reads.tar.gz
    tar zxvf reads.tar.gz

Now, the he raw reads should be quality trimmed (e.g. using sickle).

.. code-block:: sh

    sickle pe -f reads1.fastq -r reads2.fastq -t sanger -o \
    reads1.trim.fastq -p reads2.trim.fastq -s reads.singles.fastq -q 20

Given the species of interest (e.g. `P. acnes`), download and untar the 
precomputed Bowtie reference database available at 
ftp://ftp.fmach.it/metagenomics/strainest/ref/ (e.g. ``pacnes.tar.gz``):

.. code-block:: sh

    wget ftp://ftp.fmach.it/metagenomics/strainest/ref/pacnes.tar.gz
    tar zxvf pacnes.tar.gz

The Bowtie2 database is available in the ``P_acnes/bowtie`` directory. Now we can
align the metagenome against the database:

.. code-block:: sh

    bowtie2 --very-fast --no-unal -x P_acnes/bowtie/align -1 reads1.trim.fastq \
    -2 reads2.trim.fastq -S reads.sam

Now we can sort and index the BAM file:

.. code-block:: sh

    samtools view -b reads.sam > reads.bam
    samtools sort reads.bam -o reads.sorted.bam
    samtools index reads.sorted.bam

Finally, the ``strainest est`` command predicts the strain abundances:

.. code-block:: sh

    strainest est snp_clust.dgrp reads.sorted.bam outputdir

In the output directory we can find:

    abund.txt
        The predicted abundances for each reference genome

    max_ident.txt
        For each reference genome, the percentage of alleles that are present in
        the metagenome

    info.txt
        Information about the prediction, including the prediction Pearson R

    counts.txt
        Number of counts for each SNV position/base pairs

    mse.pdf
        Lasso cross-validation plot as a function of the shrinkage coefficient


(Optional) Build a custom reference SNV profile
-----------------------------------------------
See the Online Methods in the paper.