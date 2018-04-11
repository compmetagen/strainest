## This code is written by Davide Albanese, <davide.albanese@gmail.com>
## Copyright (C) 2015 Fondazione Edmund Mach
## Copyright (C) 2015 Davide Albanese

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

import warnings
import sys
import csv
import os.path
import collections

import numpy as np
import scipy as sp
from scipy import sparse
import pandas as pd
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    from sklearn.cross_validation import ShuffleSplit
from sklearn.linear_model import Lasso, LassoCV
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import pysam


def get_counts(samfile, positions, quality_threshold=15):
    """Get the pileup column at positions pos (1-based)."""

    counts = np.zeros((4, len(positions)), dtype=np.int64)
    for i, position in enumerate(positions):
        for reference in samfile.references:
            ref_count = samfile.count_coverage(
                contig=reference, start=position-1, end=position,
                quality_threshold=quality_threshold)
            counts[:, i] += np.asarray(ref_count, dtype=np.int64).flatten()
    return counts


def lasso_mpm(alphas, mse_path):
    """ Compute the Lasso most parsimonious model (i.e. fewest genomes, smaller
    alpha) within one standard errors of the best model.
    """

    mse_mean = np.mean(mse_path, axis=1)
    mse_std = np.std(mse_path, axis=1)
    mse_min_idx = np.argmin(mse_mean)
    mse_min = mse_mean[mse_min_idx]
    mse_min_std = mse_std[mse_min_idx]
    mse_min_std_min = mse_min - mse_min_std
    mse_min_std_max = mse_min + mse_min_std

    mse_mpm_idx = mse_min_idx
    for i in range(mse_min_idx-1, -1, -1):
        if (mse_mean[i]>=mse_min_std_min) and (mse_mean[i]<=mse_min_std_max):
            mse_mpm_idx = i

    alpha_mpm = alphas[mse_mpm_idx]
    mse_mean_mpm = mse_mean[mse_mpm_idx]
    mse_std_mpm = mse_std[mse_mpm_idx]

    return alpha_mpm, mse_mean_mpm, mse_std_mpm


def est(snp_fn, bam_fn, output_dir, quality_thr=20, min_depth_percentile=10,
        max_depth_percentile=90, min_depth_absolute=6, min_depth_base=0.01,
        max_ident_thr=0.95, threads=1):

    def write_abund_info():
        abund_df.to_csv(abund_fn, float_format="%.6f", sep='\t',
                        header=[bam_bn], index_label="OTU")

        with open(info_fn, 'w') as info_handle:
            info_writer = csv.writer(info_handle, delimiter='\t',
                                     lineterminator='\n')
            info_writer.writerow(info.keys())
            info_writer.writerow(info.values())


    CV_NITER = 20
    NALPHA = 50
    MAX_NITER = 5000
    TEST_SIZE = 0.5

    BASE_ORDER = ['A', 'C', 'G', 'T']

    BASE_P = {
        'A': [1, 0, 0, 0],
        'C': [0, 1, 0, 0],
        'G': [0, 0, 1, 0],
        'T': [0, 0, 0, 1],
    }

    # output filenames
    abund_fn = os.path.join(output_dir, "abund.txt")
    max_ident_fn = os.path.join(output_dir, "max_ident.txt")
    info_fn = os.path.join(output_dir, "info.txt")
    counts_fn = os.path.join(output_dir, "counts.txt")
    mse_fn = os.path.join(output_dir, "mse.pdf")

    snp = pd.read_csv(snp_fn, sep=',', index_col=0)

    snp = snp.drop('Ref', 1)
    ref_id = snp.index.name
    positions = snp.index.values # positions (1-based)
    genomes = snp.columns.values

    abund_df = pd.Series(index=genomes)
    abund_df.fillna(0.0, inplace=True)
    info = collections.OrderedDict()
    info["Filename"] = ""
    info["MinDepth"] = np.nan
    info["MaxDepth"] = np.nan
    info["NPos"] = np.nan
    info["NGen"] = np.nan
    info["Alpha"] = np.nan
    info["MSEAve"] = np.nan
    info["MSEStd"] = np.nan
    info["NLassoIter"] = np.nan
    info["R"] = np.nan
    info["PVal"] = np.nan

    # get counts
    bam_bn = os.path.basename(bam_fn)
    align = pysam.AlignmentFile(bam_fn, 'rb')
    counts = get_counts(align, positions, quality_thr)
    align.close()

    info["Filename"] = bam_bn

    # compute depths of coverage (DoC) and percentiles
    depths = counts.sum(axis=0)
    depths_nz = depths[depths>0]
    min_depth, max_depth = np.percentile(
        depths_nz, [min_depth_percentile, max_depth_percentile])
    info["MinDepth"] = min_depth
    info["MaxDepth"] = max_depth

    # filter by DoC
    keep = np.logical_and.reduce(
        (depths>=min_depth, depths<=max_depth, depths>=min_depth_absolute))
    snp = snp[keep]
    counts = counts[:, keep]
    positions_used = snp.index.values

    info["NPos"] = snp.shape[0]
    if info["NPos"] < 2:
        sys.stderr.write("Insufficient number of positions covered (<2)")
        write_abund_info()
        return

    # count filter
    counts[counts<(counts.sum(axis=0)*min_depth_base)] = 0

    # compute frequencies and alleles
    freqs = (counts / counts.sum(axis=0))
    alleles_count = (counts>0).sum(axis=0)

    # write counts
    counts_df = pd.DataFrame(index=positions_used, columns=BASE_ORDER,
                             data=counts.T)
    counts_df.to_csv(counts_fn, sep='\t', index_label="Pos")

    # SNP probability matrix
    snp_prob_flat = np.empty((4*snp.shape[0], snp.shape[1]), dtype=np.float64)
    for j in range(snp.shape[1]):
        p = [BASE_P.get(b, [0, 0, 0, 0]) for b in snp.iloc[:, j]]
        snp_prob_flat[:, j] = np.asarray(p).flatten()
    freqs_flat = freqs.flatten(order='F')

    # max identities
    tmp = snp_prob_flat * freqs_flat.reshape(-1, 1)
    tmp = tmp.reshape(-1 , 4, tmp.shape[1]).sum(axis=1)
    max_ident = (tmp>0).sum(axis=0) / tmp.shape[0]

    # write maximum identities
    max_ident_df = pd.Series(index=genomes, data=max_ident)
    max_ident_df.to_csv(max_ident_fn, float_format="%.6f", sep='\t',
                        header=[bam_bn], index_label="OTU")

    # remove genomes with the identity below the threshold
    keep = (max_ident >= max_ident_thr)
    snp_prob_flat = snp_prob_flat[:, keep]
    genomes_used = genomes[keep]
    info["NGen"] = genomes_used.shape[0]

    if info["NGen"] == 0:
        sys.stderr.write("No compatible genomes available")
        write_abund_info()
        return

    # remove position/base pairs where the sum of probabilities is 0
    keep = (snp_prob_flat.sum(axis=1) != 0)
    snp_prob_flat = snp_prob_flat[keep]
    freqs_flat = freqs_flat[keep]

    # Lasso
    X = sparse.csr_matrix(snp_prob_flat)
    y = freqs_flat

    # Tuning Lasso alpha (most parsimonious model)
    cv = ShuffleSplit(X.shape[0], n_iter=CV_NITER, test_size=TEST_SIZE,
                      random_state=0)

    lasso_cv = LassoCV(eps=0.001, n_alphas=NALPHA,
                       fit_intercept=False, normalize=False,
                       precompute='auto', max_iter=MAX_NITER,
                       tol=0.0001, copy_X=True, cv=cv, verbose=False,
                       n_jobs=threads, positive=True, random_state=0,
                       selection='cyclic')
    lasso_cv.fit(X, y)
    alpha, mse_ave, mse_std = lasso_mpm(lasso_cv.alphas_, lasso_cv.mse_path_)

    info["Alpha"] = alpha
    info["MSEAve"] = mse_ave
    info["MSEStd"] = mse_std

    # MSE plot
    m_log_alphas = -np.log10(lasso_cv.alphas_)
    fig = plt.figure(1)
    plt.plot(m_log_alphas, lasso_cv.mse_path_, ':')
    plt.plot(m_log_alphas, lasso_cv.mse_path_.mean(axis=-1), 'k',
             label='Average across the folds', linewidth=2)
    plt.axvline(-np.log10(alpha), linestyle='--', color='k',
                label='Alpha (most parsimonious model)')
    plt.xlabel('-log(alpha)')
    plt.ylabel('Mean square error')
    fig.savefig(mse_fn, bbox_inches='tight')

    # Relative abundance estimation
    lasso = Lasso(alpha=alpha, fit_intercept=False, normalize=False,
                  precompute=False, copy_X=True, max_iter=MAX_NITER,
                  tol=0.0001, warm_start=False, positive=True,
                  random_state=0, selection='cyclic')
    lasso.fit(X, y)
    lasso_coef = np.atleast_1d(lasso.coef_)

    info["NLassoIter"] = lasso.n_iter_

    # normalize the Lasso coefficients
    coef_norm = lasso_coef / np.sum(lasso_coef)
    abund_df[genomes_used] = coef_norm

    # pearson correlation
    y_pred = X.dot(coef_norm)

    info["R"], info["PVal"] = pearsonr(y, y_pred)

    write_abund_info()

    if (np.max(coef_norm) < 0.9) and (max_ident_thr < 0.99):
        sys.stdout.write("WARNING: the maximum identity threshold is <0.99 and "
                         "StrainEst has inferred a mixture of strains. The "
                         "mixture of strains could be a single strain with no "
                         "available reference genome. Please check the file "
                         "counts.txt.")
