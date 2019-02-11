#!/usr/bin/env python
# -*- coding: utf-8

import sys
import logging

from wiz.qc import gc


def run(args):
    """
    main function for wiz quality check (qc)
    """
    logger = logging.getLogger(__name__)

    try:
        gc_per_bin = []
        for genome in args.genome:
            # TODO parse with Bio.SeqIO
            average_gc = gc.average_gc(seq, truncate=False)
            gc_per_bin.append((average_gc, name))

            gc.scatter_gc(gc_per_bin)
            gc.distplot_gc(gc_per_bin)
            gc.histogram_gc(gc_per_bin)
    except expression as identifier:
        logger.error("something bad happened")
        sys.exit(1)
    else:
        logger.info("wiz qc finished. Goodbye.")
