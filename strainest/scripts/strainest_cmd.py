#! /usr/bin/env python

##    Copyright 2016 Davide Albanese <davide.albanese@gmail.com>
##    Copyright 2016 Fondazione Edmund Mach (FEM)

##    This file is part of strainest.
##
##    strainest is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    strainest is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.

##    You should have received a copy of the GNU General Public License
##    along with strainest.  If not, see <http://www.gnu.org/licenses/>.

import os
import os.path

import click
from Bio import SeqIO

import strainest.api
from strainest import __version__


@click.group()
def cli():
    """StrainEst - abundance estimation of strains
    """
    pass

@cli.command()
@click.argument('genomes', required=True, nargs=-1,
                type=click.Path(exists=True))
@click.argument('reference', type=click.Path(exists=True))
@click.argument('mapped', type=click.Path(exists=False, writable=True))
def mapgenomes(genomes, reference, mapped):
    """Map genomes to a reference.

    Align one or more genomes to a reference genome. Only the first sequence
    in the reference genome is considered. Input and output files must be in
    FASTA format.

    EXAMPLES

        strainest mapgenome genome1.fna genome2.fna reference.fna mapped.fna
    """
    with open(mapped, 'wb') as mapped_handle:
        with click.progressbar(genomes) as bar:
            for genome in bar:
                record = strainest.api.mapgenome(genome, reference, mapped)
                SeqIO.write(record, mapped_handle, "fasta")


@cli.command()
@click.argument('mapped', type=click.Path(exists=True))
@click.argument('output', type=click.Path(writable=True, dir_okay=True,
                resolve_path=True))
def mapstats(mapped, output):
    """Compute basic statistics on the mapped genomes file.

    EXAMPLES

        strainest mapstats mapped.fna stats
    """
    try:
        os.makedirs(output)
    except OSError:
        if not os.path.isdir(output):
            click.echo("directory {} cannot be created".format(output),
                       err=True)
            exit(1)

    strainest.api.mapstats(mapped, output)


@cli.command()
@click.argument('reference', type=click.Path(exists=True))
@click.argument('mapped', type=click.Path(exists=True))
@click.argument('output', type=click.Path(exists=False, writable=True))
def map2snp(reference, mapped, output):
    """Build the SNP matrix in DGRP format.

    EXAMPLES

        strainest map2snp reference.fna mapped.fna snp.dgrp
    """

    strainest.api.map2snp(reference, mapped, output)

@cli.command()
@click.argument('snp', type=click.Path(exists=True))
@click.argument('dist', type=click.Path(exists=False, writable=True))
@click.argument('hist', type=click.Path(exists=False, writable=True))
def snpdist(snp, dist, hist):
    """Compute the Hamming distances between sequences in SNP matrix (in DGRP
    format). Moreover, it returns the distances histogram in HIST.

    EXAMPLES

        strainest snpdist snp.dgrp dist.txt hist.pdf
    """

    strainest.api.snpdist(snp, dist, hist)


@cli.command()
@click.argument('snp', type=click.Path(exists=True))
@click.argument('dist', type=click.Path(exists=True))
@click.argument('snpout', type=click.Path(exists=False, writable=True))
@click.argument('clust', type=click.Path(exists=False, writable=True))
@click.option('-t', '--thr', type=click.FLOAT, default=0.01,
              show_default=True, help="distance threshold")
def snpclust(snp, dist, snpout, clust, thr):
    """Given a SNP matrix (in DGRP format) and a distance matrix,
    snpdist clusters the profiles returning a SNP matrix containing the
    representative profiles (SNPOUT) and a cluster file (CLUST).

    EXAMPLES
        strainest map2snp reference.fna mapped.fna snp.dgrp
    """

    strainest.api.snpclust(snp, dist, snpout, clust, thr)


@cli.command()
@click.argument('snp', type=click.Path(exists=True))
@click.argument('bam', type=click.Path(exists=True))
@click.argument('output', type=click.Path(writable=True, dir_okay=True,
                resolve_path=True))
@click.option('-q', '--quality-thr', type=click.INT, default=20,
              show_default=True, help="base quality threshold")
@click.option('-p', '--min-depth-percentile', type=click.FLOAT, default=10,
              show_default=True, help="discard positions where the depth of "
              "coverage is lower than the MIN_DEPTH_PERCENTILE percentile")
@click.option('-P', '--max-depth-percentile', type=click.FLOAT, default=90,
              show_default=True, help="discard positions where the depth of "
              "coverage is higher than the MAX_DEPTH_PERCENTILE percentile")
@click.option('-a', '--min-depth-absolute', type=click.INT, default=6,
              show_default=True, help="discard positions where the depth of "
              "coverage is lower than the MIN_DEPTH_ABSOLUTE")
@click.option('-b', '--min-depth-base', type=click.FLOAT, default=0.01,
              show_default=True, help="filter base counts (set to 0) where "
              "they are lower then MIN_DEPTH_BASE x DoC (applied independently 
              "for each allelic position)")
@click.option('-d', '--max-ident-thr', type=click.FLOAT, default=0.95,
              show_default=True, help="discard genomes with less than "
              "MAX_IDENT_THR maximum identity")
@click.option('-t', '--threads', type=click.INT, default=1,
              show_default=True, help="number of threads to use in model "
              "selection")
def est(snp, bam, output, quality_thr, min_depth_percentile,
        max_depth_percentile, min_depth_absolute, min_depth_base, max_ident_thr,
        threads):
    """Estimate relative abundance of strains. The BAM file must be sorted and
    indexed.

    EXAMPLES

        strainest est snp.dgrp align.bam align -t 4
    """


    try:
        os.makedirs(output)
    except OSError:
        if not os.path.isdir(output):
            click.echo("directory {} cannot be created".format(output),
                       err=True)
            exit(1)

    strainest.api.est(snp, bam, output, quality_thr,
                      min_depth_percentile, max_depth_percentile,
                      min_depth_absolute, min_depth_base, max_ident_thr,
                      threads)
