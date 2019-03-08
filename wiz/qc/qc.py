#!/usr/bin/env python
# -*- coding: utf-8

import sys
import logging

from wiz.qc import report
from wiz.qc.bins import Bins
from wiz.annotate.tools import prodigal
from wiz.qc.tools import get_file_name as gfn
from wiz.qc.tools import get_gene_seq as ggs
from Bio import SeqIO
import os


def run(args):
    """
    main function for wiz quality check (qc)
    """
    logger = logging.getLogger(__name__)
    logger.debug(" [START] wiz qc")

    try:
        bins = []
        report.create_dir(args.output)
        for genome_file in args.genomes:
            logger.debug(f"Prodigal on {genome_file}")
            prodigal(genome_file,output_dir=args.output,output_prefix=gfn(genome_file))
            seq_genes = ggs(args.output,gfn(genome_file))
            logger.debug(f"genome file: {genome_file}")
            f = open(genome_file, 'r')
            with f:
                fasta_file = SeqIO.parse(f, 'fasta')
                logger.info(f"Loading {genome_file}")
                for record in fasta_file:
                    genome_bin = Bins(record, args.window, seq_genes)
                    bins.append(genome_bin)

        report_data = report.Report(bins, args.window)
        report_html = report.jinja_report(report_data, args)
        report.write_QCreport(args, report_html)
    except ValueError as Ve:
        logger.error(" something bad happened")
        logger.error(Ve)
        sys.exit(1)
    except KeyboardInterrupt:
        logger.error(" something bad happened ?")
        logger.error(" You cancelled the operations.")
    else:
        logger.info(" wiz qc finished. Goodbye.")

# TODO adding parameter
# TODO converting genome parameter in genome folder
