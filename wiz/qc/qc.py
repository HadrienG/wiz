#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import logging
from Bio import SeqIO
from wiz.qc import report, tools, taxonomy
from wiz.qc.bins import Bins, Contig
from wiz.qc.phylo import phylogeny
from wiz.annotate.tools import prodigal


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
            f = open(genome_file, 'r')
            with f:
                fasta_file = SeqIO.parse(f, 'fasta')
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
                if contig.id in genome_genes.keys():
                    contig.genes_position = genome_genes[contig.id]
                    contig.coding_density = taxonomy.coding_density(
                        contig.genes_position,
                        contig.length)
                # if the sequence in fasta file don't have a correct name
                # for prodigal, we try with the var "second_name" who contain
                # a name like prodigal default name
                elif second_name in genome_genes.keys():
                    contig.genes_position = genome_genes[second_name]
                    contig.coding_density = taxonomy.coding_density(
                        contig.genes_position,
                        contig.length)
            # phylogenic analysis based on prodigal output
            phylo = phylogeny(b.filename, args)
            for contig in b.contigs:
                if contig.id in phylo.keys():
                    contig.phylogeny_profils = phylo[contig.id]
            # taxonomic analysis
            # for contig in b.contigs:
            #     contig.taxonomy = taxonomy.taxonomy(
            #         contig.id,
            #         contig.sequence,
            #         args)

        report_data = report.Report(bins, args)
        report.jinja_report(report_data, args)

    except ValueError as Ve:
        logger.error(" something bad happened")
        logger.error(Ve)
        sys.exit(1)
    except KeyboardInterrupt:
        logger.error(" something bad happened ?")
        logger.error(" You cancelled the operations.")
    else:
        logger.info(" wiz qc finished. Goodbye.")

# TODO search this line in each file
# TODO " # TODO uncomment before release "
