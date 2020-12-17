#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#  Copyright 2019 Luca Beltrame <luca.beltrame@marionegri.it>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with utils. If not, see <http://www.gnu.org/licenses/>.

from contextlib import contextmanager
import tempfile
import os
from pathlib import Path

from arghandler import subcmd, ArgumentHandler

from joblib import delayed, Parallel
from natsort import natsorted
import pysam

CHRS = ["chr1", "chr2", "chr3", "chr4",
        "chr5", "chr6", "chr7", "chr8",
        "chr9", "chr10", "chr11", "chr12",
        "chr13", "chr14", "chr15", "chr16",
        "chr17", "chr18", "chr19", "chr20",
        "chr21", "chr22", "chrX"]


@contextmanager
def cd(subpath):
    old_path = Path.cwd()
    os.chdir(str(old_path / subpath))
    try:
        yield
    finally:
        os.chdir(str(old_path))


def size_select_bam(source_file, chromosome, destination: Path, cores=4,
                    index_filename=None):

    srcbam = pysam.AlignmentFile(str(source_file), threads=cores/2,
                                 index_filename=str(index_filename))
    destination = (destination / (chromosome+"_"+source_file.name))

    with pysam.AlignmentFile(str(destination), template=srcbam, mode="wb",
                             threads=cores/2) as handle:
        for read in srcbam.fetch(chromosome, multiple_iterators=False):

            if read.qlen < 90 or read.qlen > 150:
                continue

            handle.write(read)

    return str(destination)


def parallel_computation(source_file: Path, destination: Path, cores=4,
                         index_source=False):

    with tempfile.TemporaryDirectory() as temp:
        temp2 = Path(temp)
        processor = delayed(size_select_bam)
        pool = Parallel(n_jobs=cores)

        index_path = Path.cwd() / (source_file.name + ".bai")
        pysam.index("-@", str(cores), str(source_file), str(index_path))
        assert index_path.exists()
        result = pool(processor(source_file, chrom, temp2, cores, index_path)
                      for chrom in CHRS)
        result = natsorted(result, key=lambda x: x.split(".")[0])
        pysam.merge("-h", str(source_file), "-@", str(cores), str(destination),
                    *result)
        pysam.index("-@", str(cores), str(destination))


@subcmd
def single(parser, context, args):

    parser.add_argument("-c", "--cores", type=int,
                        help="Number of cores to use")
    parser.add_argument(
        "--index-source",
        action="store_true",
        help="Index source BAM in current directory (use with Nextflow)"
    )
    parser.add_argument("source", help="Source BAM file")
    parser.add_argument("destination", help="Destination BAM to create")
    options = parser.parse_args(args)
    source = Path(options.source).absolute()
    destination = Path(options.destination).absolute()

    result = parallel_computation(source, destination, cores=options.cores,
                                  index_source=options.index_source)


def main():

    handler = ArgumentHandler()
    handler.run()


if __name__ == '__main__':
    main()
