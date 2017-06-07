## This code is written by Davide Albanese, <davide.albanese@gmail.com>
## Copyright (C) 2015-2016 Fondazione Edmund Mach
## Copyright (C) 2015-2016 Davide Albanese

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division

import os
import os.path

import matplotlib as mpl
mpl.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser


def mapstats(mapped_fn, output_dir):

    CORE_FN = os.path.join(output_dir, "core.txt")
    COVERAGE_FN = os.path.join(output_dir, "coverage.txt")
    CORE_PLOT_FN = os.path.join(output_dir, "core_plot.pdf")
    COVERAGE_PLOT_FN = os.path.join(output_dir, "coverage_plot.pdf")

    # read the sequence length
    with open(mapped_fn, 'rb') as mapped_handle:
        title, seq = SimpleFastaParser(mapped_handle).next()
        seqlen = len(seq)

    # load the number of bases mapping for each sequence
    mapinfo = []
    with open(mapped_fn, 'rb') as mapped_handle:
        for title, seq in SimpleFastaParser(mapped_handle):
            if len(seq) != seqlen:
                raise ValueError("Sequences in the mapped file must have the "
                                 "same length")
            seqid = title.split()[0]
            mapn = seqlen - seq.count('-')
            mapinfo.append([seqid, mapn])

    # sort mapinfo by the number of bases mapping in decreasing order
    mapinfo_sorted = sorted(mapinfo, key=lambda elem: elem[1], reverse=True)
    seqn = len(mapinfo)

    map_records = SeqIO.index(mapped_fn, 'fasta')
    coverage, core_list = np.zeros(seqlen, dtype=np.int), []
    for i, (seqid, mapn) in enumerate(mapinfo_sorted):
        coverage += (np.asarray(list(map_records[seqid].seq)) != '-')
        core_list.append([seqid, mapn, (coverage == i+1).sum()])

    # write coverage and core files
    coverage_df = pd.DataFrame(
        coverage, index=range(1, coverage.shape[0]+1), columns=["Cov"])
    coverage_df.index.name = "Pos"
    coverage_df.to_csv(COVERAGE_FN, sep='\t', header=True)

    core_df = pd.DataFrame(core_list, columns=["SeqID", "MapN", "CoreN"])
    core_df.set_index("SeqID", inplace=True)
    core_df.to_csv(CORE_FN, sep='\t', header=True)

    # plot coverage
    fig = plt.figure(figsize=(12, 6))
    ax = plt.subplot(111)
    coverage_count = np.bincount(coverage_df["Cov"], minlength=seqn+1)[1:]
    plt.bar(np.arange(1, seqn+1)-0.5, coverage_count, width=1, log=False,
            linewidth=0, color="#0072B2")
    plt.title("Coverage histogram")
    plt.xlabel("Coverage")
    plt.ylabel("# of bases")
    plt.xlim((0.5, seqn+0.5))
    plt.ylim((1, ax.get_ylim()[1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.savefig(COVERAGE_PLOT_FN, bbox_inches='tight')

    # plot # bases mapped vs # of bases covered (core)
    fig = plt.figure(figsize=(8, 8))
    ax = plt.subplot(111)
    plt.plot(core_df["MapN"], core_df["CoreN"], marker='o',
             markerfacecolor="#0072B2", markeredgewidth=0,
             markersize=8, linewidth=2, color="#0072B2")
    plt.xlabel("# of bases mapping")
    plt.ylabel("# of bases core")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.savefig(CORE_PLOT_FN, bbox_inches='tight')
