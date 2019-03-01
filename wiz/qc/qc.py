#!/usr/bin/env python
# -*- coding: utf-8

import sys
import logging

from wiz.qc import report
from wiz.qc.bin import Bin
from Bio import SeqIO


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
                    genome_bin = Bin(record, args.window)
                    bins.append(genome_bin)

        report_data = report.Report(bins, args.window)
        report_html = report.jinja_report(report_data)
        with open(args.output, "w") as r:
            r.writelines(report_html)
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
