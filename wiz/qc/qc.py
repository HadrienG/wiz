#!/usr/bin/env python
# -*- coding: utf-8

import sys
import logging

from wiz.qc import gc, tools, bins, tetra, report


def run(args):
    """
    main function for wiz quality check (qc)
    """
    logger = logging.getLogger(__name__)
    logger.debug("[START] wiz qc")

    try:
        genomes = []
        logger.debug("Start import")
        for gen in args.genome:
            logger.debug(f"genome file: {gen}")
            genomes.append(tools.genome_parser(gen))
        logger.debug("Import ok")
        for gen in genomes:
            bin = bins.Bins(gen)
            bin.gc = gc.average_gc(bin.seq, truncate=True)
            bin.tetra = tetra.tetranuc_count(bin.seq)
            bin.gc = gc.percentil_filter(bin.gc)
            report.scatter_gc(bin.gc)
            report.distplot_gc(bin.gc)
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
