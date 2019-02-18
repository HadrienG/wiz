#!/usr/bin/env python
# -*- coding: utf-8

import sys
import logging

from wiz.qc import gc, tools, bins, tetra


def run(args):
    """
    main function for wiz quality check (qc)
    """
    logger = logging.getLogger(__name__)
    logger.info("wiz qc started !")

    try:
        gc_per_bin = []
        logger.info(f"genome: {args.genome}")
        bin = bins.Bins(tools.genome_parser(args.genome))
        bin.gc = gc.average_gc(bin.seq, truncate=True)
        bin.tetra = tetra.tetranuc_count(bin.seq)
        bin.gc = gc.percentil_filter(bin.gc)
        gc.scatter_gc(bin.gc)
        gc.distplot_gc(bin.gc)
    except ValueError as Ve:
        logger.error("something bad happened")
        logger.error(Ve)
        sys.exit(1)
    except KeyboardInterrupt:
        logger.error("something bad happened ?")
        logger.error("You cancelled the operation.")
    else:
        logger.info("wiz qc finished. Goodbye.")

# TODO adding parameter
# TODO converting genome parameter in genome folder
