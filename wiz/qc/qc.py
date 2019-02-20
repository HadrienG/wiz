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
        for gen in args.genome:
            logger.debug(f"genome file: {gen}")
            genomes += (tools.genome_parser(gen, auto_import_folder=args.folder))
        bins_list = []
        if len(genomes) == 0:
            logger.error("No sequence provided")
            exit(2)
        else:
            for gen in genomes:
                logger.debug(f"gen : {gen}")
                bin = bins.Bins(gen)
                bins_list.append(bin)
            del genomes
            if not args.ignorefilter:
                bins_list = tools.check_for_duplicates(bins_list, args.autofilter)
            for bin in bins_list:
                bin.gc = gc.average_gc(bin.seq, window_size=args.window, truncate=True)
                bin.gc_filtered = gc.percentil_filter(bin.gc, bin.gc_percentil)
                bin.gc_bounds = gc.get_bounds(bin.gc, bin.gc_percentil)
                bin.tetra = tetra.tetranuc_count(bin.seq)
            report_data = report.Report(bins_list, args.window)
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
