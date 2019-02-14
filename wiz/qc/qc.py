#!/usr/bin/env python
# -*- coding: utf-8

import sys
import logging

from wiz.qc import gc, tools, qcclass, tetra


def run(args):
    """
    main function for wiz quality check (qc)
    """
    logger = logging.getLogger(__name__)

    try:
        gc_per_bin = []
        for genome in args.genome:
            seq = qcclass.bins(tools.genome_parser(genome, "fasta"))
            seq.gc = gc.average_gc(seq.seq, truncate=True)
            seq.tetra = tetra.tetranuc_count(seq.seq)
            gc.scatter_gc(seq.gc)
            gc.distplot_gc(seq.gc)
            gc.histogram_gc(seq.gc)
    except expression as identifier:
        logger.error("something bad happened")
        sys.exit(1)
    else:
        logger.info("wiz qc finished. Goodbye.")

#TODO adding parameter
