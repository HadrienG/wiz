#!/usr/bin/env python
# -*- coding: utf-8

"""
main function of module QC.
"""

import sys
import logging

from Bio import SeqIO

from wiz.qc import report
from wiz.qc import tools
# from wiz.qc import taxonomy
from wiz.qc import cd
from wiz.qc.bins import Bins
from wiz.qc.bins import Contig
from wiz.annotate.tools import prodigal
# from wiz.qc.phylo import phylogeny


def run(args):
    """
    main function for wiz quality check (qc)
    """
    logger = logging.getLogger(__name__)
    logger.debug(" Welcome in wiz qc")
    try:
        bins = []
        tools.create_dirs(args.output, args.force)

        for genome_file in args.genomes:
            file = open(genome_file, 'r')
            with file:
                fasta_file = SeqIO.parse(file, 'fasta')
                logger.info(f" Loading {genome_file}")
                contigs = []

                for record in fasta_file:
                    contigs.append(Contig(record, args.window))

                genome_bin = Bins(genome_file, contigs)
                bins.append(genome_bin)

        for b in bins:
            # coding density analysis
            logger.debug(f" Run prodigal for {b.filename}")
            output, prefix = args.output+"/prodigal", b.filename
            prodigal(b.path, output_dir=output, output_prefix=prefix)
            genome_genes = tools.genes(output, b.filename)

            for contig in b.contigs:
                second_name = f"Prodigal_Seq_{str(b.contigs.index(contig)+1)}"
                if contig.uid in genome_genes.keys():
                    contig.genes_position = genome_genes[contig.uid]
                    contig.coding_density = cd.coding_density(
                        contig.genes_position,
                        contig.length)

                # if the sequence in fasta file don't have a correct name
                # for prodigal, we try with the var "second_name" who contain
                # a name like prodigal default name
                elif second_name in genome_genes.keys():
                    contig.genes_position = genome_genes[second_name]
                    contig.coding_density = cd.coding_density(
                        contig.genes_position,
                        contig.length)

            # # taxonomic analysis
            # TODO uncomment before release
            # for contig in b.contigs:
            #     contig.taxonomy = taxonomy.taxonomy(
            #         contig.uid,
            #         contig.sequence,
            #         args)

            # # Implementation of phylogenic analysis in the bad module
            # # phylogenic analysis based on prodigal output
            # phylo = phylogeny(b.filename, args)
            # for contig in b.contigs:
            #     if contig.uid in phylo.keys():
            #         contig.phylogeny_profils = phylo[contig.uid]

        report_data = report.Report(bins, args)
        report.jinja_report(report_data, args)

    except ValueError as error:
        logger.error(" something bad happened")
        logger.error(error)
        sys.exit(1)
    except KeyboardInterrupt:
        logger.error(" something bad happened ?")
        logger.error(" You cancelled the operations.")
    else:
        logger.info(" wiz qc finished. Goodbye.")
