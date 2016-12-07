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

import re


class Alignment:
    pass

class MummerAlignmentReader:
    def __init__(self, handle):
        self.__handle = handle

        # load filenames
        line = self.__handle.readline()
        filenames = line.split()
        self.filename1 = filenames[0]
        self.filename2 = filenames[1]

        # load sequence ids
        while True:
            line = self.__handle.readline()
            if line.startswith('-- Alignments between'):
                ids = line.split()
                self.id1 = ids[3]
                self.id2 = ids[5]
                break

    def __iter__(self):
        return self

    def next(self):

        # read alignment header
        while True:
            line = self.__handle.readline()
            if not line:
                raise StopIteration
            if line.startswith('-- BEGIN alignment'):
                coords = re.findall('[+-]?\d+', line)
                coords = [int(elem) for elem in coords]
                alignment = Alignment()
                alignment.strand1 = '+' if coords[0] == +1 else '-'
                alignment.start1 = coords[1]
                alignment.end1 = coords[2]
                alignment.strand2 = '+' if coords[3] == +1 else '-'
                alignment.start2 = coords[4]
                alignment.end2 = coords[5]
                break

        # read sequences
        seq1, seq2 = '', ''
        while True:
            line = self.__handle.readline()
            mseq1 = re.search(r'\d+ +[a-zA-Z.]+', line)
            if mseq1:
                seq1 += mseq1.group().split()[1]
                line = self.__handle.readline()
                seq2 += re.search(r'[a-zA-Z.]+', line).group()
            if line.startswith('--   END alignment'):
                alignment.seq1 = seq1
                alignment.seq2 = seq2
                break

        return alignment
