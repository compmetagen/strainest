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

import argparse
import csv
import textwrap

import numpy as np
from scipy.cluster import hierarchy
import pandas as pd
from scipy.spatial.distance import squareform


def snpclust(snp_fn, dist_fn, snpout_fn, clust_fn, thr=0.01):

    # read SNP and distance matrix files
    snp = pd.read_csv(snp_fn, sep=',', index_col=0)
    dist = pd.read_csv(dist_fn, sep=',', index_col=0)

    # hierarchical clustering
    Z = hierarchy.complete(squareform(dist))
    clust_ids = hierarchy.fcluster(Z, t=thr, criterion='distance')

    # compute the cluster representatives and build the cluster dictionary
    clust = dict()
    for ci in np.unique(clust_ids):
        idx = np.where(ci==clust_ids)[0]
        if idx.shape[0] == 1:
            r = dist.index[idx[0]]
            clust[r] = [r]
        else:
            dist_sum = dist.iloc[idx, idx].sum(axis=0)
            clust[dist_sum.idxmin()] = dist.index[idx].tolist()

    # write the SNP output file containing the medoids only
    snp_out = pd.concat((snp["Ref"], snp[clust.keys()]), axis=1)
    snp_out.to_csv(snpout_fn, index_label=snp.index.name)

    # write the cluster file
    with open(clust_fn, 'wb') as clust_handle:
        clust_writer = csv.writer(clust_handle, delimiter='\t',
                                  lineterminator='\n')
        for r, ms in clust.items():
            for m in ms:
                clust_writer.writerow([m, r, "{:.6f}".format(dist.loc[m, r])])
