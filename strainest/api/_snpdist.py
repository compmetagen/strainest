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

import click
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform


def snpdist(snp_fn, dist_fn, hist_fn):

    snp = pd.read_csv(snp_fn, sep=',', index_col=0)
    snp.drop('Ref', axis=1, inplace=True)
    dist = pd.DataFrame(index=snp.columns, columns=snp.columns, dtype=np.float)
    dist.fillna(0.0, inplace=True)
    distn = int((snp.shape[1] * (snp.shape[1]-1)) / 2)
    with click.progressbar(length=distn) as bar:
        for i in range(snp.shape[1]):
            for j in range(i+1, snp.shape[1]):
                d = np.sum(snp.iloc[:, i] != snp.iloc[:, j]) / snp.shape[0]
                dist.iloc[i, j] = d
                dist.iloc[j, i] = d
                bar.update(1)

    dist.to_csv(dist_fn, sep=',', float_format="%.6f")

    fig = plt.figure(1, figsize=(8, 4))
    plt.hist(squareform(dist), bins=100)
    plt.xlabel("Hamming distance")
    plt.ylabel("Count")
    fig.savefig(hist_fn, bbox_inches='tight')
