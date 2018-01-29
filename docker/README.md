# StrainEst - Docker

StrainEst is a reference-based method that leverages on the accumulated knowledge of the genomic variability of species of interests to identify individual strains and quantify their relative abundances in mixed metagenomic samples.

Sickle, Bowtie2 and samtools are preinstalled.

Homepage: https://github.com/compmetagen/strainest

## Available Tags/Versions

- latest: GitHub snapshot (master)
- 1.2.2: strainest 1.2.2
- 1.2.1: strainest 1.2.1


## Quickstart

1. Download the latest version:

   `docker pull compmetagen/strainest`

2. Run an instance of the image, mounting the host working directory
   (e.g. ``/Users/davide/strainest``) on to the container working directory
   ``/strainest``:

   `docker run --rm -t -i -v /Users/davide/strainest:/strainest -w /strainest compmetagen/strainest /bin/bash`

   You need to write something like ``-v //c/Users/davide/strainest:/strainest`` if
   you are in Windows or ``-v /home/davide/strainest:/strainest`` in Linux. The
   ``--rm`` option automatically removes the container when it exits.

3. Run micca without parameters:

   `root@68f6784e1101:/strainest# strainest --help`
