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

import subprocess
import os
import os.path
import re

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from strainest import mummer_path
import strainest.mummer


nucmer_path = os.path.join(mummer_path(), 'nucmer')
show_aligns_path = os.path.join(mummer_path(), 'show-aligns')


def mapgenome(genome_fn, reference_fn, delta_prefix):

    NO_ALIGN_STR = "ERROR: Could not find any alignments for"

    # run nucmer
    cmd = [nucmer_path, '--prefix', delta_prefix, reference_fn, genome_fn]
    devnull = open(os.devnull, 'w')
    proc = subprocess.Popen(cmd, stdout=devnull, stderr=subprocess.PIPE)
    _, out_stderr = proc.communicate()
    devnull.close()
    if proc.returncode:
        raise Exception(out_stderr)

    delta_fn = delta_prefix + '.delta'

    # get the first reference sequence only
    ref_record = SeqIO.parse(reference_fn, 'fasta').next()
    ref_seqid = ref_record.id
    ref_len = len(ref_record)

    genome_seq_mapped = np.asarray(['.'] * ref_len)
    genome_count_mapped = np.zeros(ref_len, dtype=np.int)
    for record in SeqIO.parse(genome_fn, 'fasta'):
        genome_seqid = record.id

        # run show-aligns
        align_fn = delta_prefix + '.align'
        with open(align_fn, 'wb') as align_handler:
            cmd = [show_aligns_path, delta_fn, ref_seqid, genome_seqid]
            proc = subprocess.Popen(cmd, stdout=align_handler,
                                    stderr=subprocess.PIPE)
            _, out_stderr = proc.communicate()

        if proc.returncode:
            if out_stderr.startswith(NO_ALIGN_STR):
                os.remove(align_fn)
                continue
            else:
                raise Exception(out_stderr)

        # load alignments
        with open(align_fn, 'rb') as align_handler:
            mar = strainest.mummer.MummerAlignmentReader(align_handler)
            for al in mar:
                ref_seq = np.asarray(list(al.seq1))
                genome_seq = np.asarray(list(al.seq2))
                genome_seq = genome_seq[np.where(ref_seq != '.')]
                genome_seq_mapped[al.start1-1:al.end1] = genome_seq
                genome_count_mapped[al.start1-1:al.end1] += 1

        os.remove(align_fn)

    os.remove(delta_fn)

    # remove repeats
    genome_seq_mapped[genome_count_mapped > 1] = '.'

    # sub '.' with '-'
    genome_seq_mapped[genome_seq_mapped == '.']  = '-'

    # write mapped genome
    out_seq = Seq(''.join(genome_seq_mapped)).upper()
    out_id = re.sub('\s+', '_', os.path.basename(genome_fn))
    out_record = SeqRecord(out_seq, id=out_id, description='')

    return out_record
