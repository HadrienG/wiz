#!/usr/bin/env python
# -*- coding: utf-8

import sys
import logging

from wiz.qc import gc, tools, bins, tetra, report
import timeit


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
                bin = bins.Bins(gen, args.window)
                bins_list.append(bin)
            # del genomes
            # if not args.ignorefilter:
            #     bins_list = tools.check_for_duplicates(bins_list, args.autofilter)
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
