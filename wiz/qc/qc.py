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
        logger.debug("[START] import")
        for gen in args.genome:
            logger.debug(f"genome file: {gen}")
            genomes += (tools.genome_parser(gen))
        logger.debug("[STOP] import")
        logger.debug("[START] Reading genome imported")
        logger.debug(genomes)
        bins_list = []
        for gen in genomes:
            logger.debug(f"gen : {gen}")
            bin = bins.Bins(gen)
            bins_list.append(bin)
        del genomes
        bins_list = tools.check_for_duplicates(bins_list)
        for bin in bins_list:
            bin.gc = gc.average_gc(bin.seq, truncate=True)
            bin.gc, bin.gc_bounds = gc.percentil_filter(bin.gc)
            bin.tetra = tetra.tetranuc_count(bin.seq)
        report.scatter_gc(bins_list)
        report.distplot_gc(bins_list)
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
