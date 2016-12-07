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

import csv

import click
import numpy as np
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser


def map2snp(reference_fn, mapped_fn, output_fn):

    # read the first sequence in reference
    reference_record = SeqIO.parse(reference_fn, 'fasta').next()
    reference_id = reference_record.id
    reference_seq = reference_record.seq.upper()
    reference_seq_arr = np.asarray(list(reference_seq))
    reference_seq_len = len(reference_seq)

    mapped_handle = open(mapped_fn, 'rb')

    seqn = 0
    coverage = np.zeros(reference_seq_len, dtype=np.int)
    click.echo("Find the core...")
    for record in SimpleFastaParser(mapped_handle):
        mapped_seq_arr = np.asarray(list(record[1]))
        coverage += (mapped_seq_arr != '-')
        seqn += 1
    core_idx = np.where(coverage == seqn)[0]
    mapped_handle.seek(0)

    snp_idx = set()
    with click.progressbar(SimpleFastaParser(mapped_handle),
                           length=seqn, label="Find SNPs") as bar:
        for record in bar:
            mapped_seq_arr = np.asarray(list(record[1]))
            diff = mapped_seq_arr[core_idx] != reference_seq_arr[core_idx]
            for idx in core_idx[diff]:
                snp_idx.add(idx)
    snp_idx = sorted(list(snp_idx))
    mapped_handle.seek(0)

    record_ids, snp = [], []
    with click.progressbar(SimpleFastaParser(mapped_handle),
                           length=seqn, label="Build the SNP matrix") as bar:
        for record in bar:
            record_ids.append(record[0].split()[0])
            mapped_seq_arr = np.asarray(list(record[1]))
            snp.append(mapped_seq_arr[snp_idx].tolist())
    snp = np.array(snp).T

    mapped_handle.close()

    click.echo("Write the snp matrix...")
    with open(output_fn, 'wb') as output_handle:
        output_writer = csv.writer(output_handle, delimiter=',',
                                   lineterminator='\n')
        output_writer.writerow([reference_id, 'Ref'] + record_ids)
        for i, s in zip(snp_idx, snp):
            output_writer.writerow([i+1, reference_seq_arr[i]] + list(s))
