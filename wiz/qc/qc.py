#!/usr/bin/env python
# -*- coding: utf-8

import sys
import logging

from Bio import SeqIO

from wiz.qc.bin import Bin
from wiz.qc import gc, tools, tetra, report


def run(args):
    """
    main function for wiz quality check (qc)
    """
    logger = logging.getLogger(__name__)
    logger.debug("[START] wiz qc")

    try:
        bins = []
        for genome_file in args.genomes:
            logger.debug(f"genome file: {genome_file}")
            f = open(genome_file, 'r')
            with f:
                fasta_file = SeqIO.parse(f, 'fasta')
                logger.info(f"Loading {genome_file}")
                for record in fasta_file:
                    genome_bin = Bin(record)
                    bins.append(genome_bin)

        for genome_bin in bins:
            genome_bin.gc = gc.average_gc(
                genome_bin.seq, window_size=args.window, truncate=True)
            genome_bin.gc_filtered = gc.percentil_filter(
                genome_bin.gc, genome_bin.gc_percentil)
            genome_bin.gc_bounds = gc.get_bounds(
                genome_bin.gc, genome_bin.gc_percentil)
            genome_bin.tetra = tetra.tetranuc_count(genome_bin.seq)
        report_data = report.Report(bins, args.window)
    except ValueError as Ve:
        logger.error("something bad happened")
        logger.error(Ve)
        sys.exit(1)
    except KeyboardInterrupt:
        logger.error("something bad happened ?")
        logger.error("You cancelled the operations.")
    else:
        logger.info("wiz qc finished. Goodbye.")

# TODO adding parameter
# TODO converting genome parameter in genome folder
