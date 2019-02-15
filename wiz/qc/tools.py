#!/usr/bin/env python
# -*- coding: utf-8

from Bio.SeqIO import parse
import os
import logging


format_supported = ["fasta"]
logger = logging.getLogger(__name__)


def genome_parser(pathfile, format="fasta"):
    if format in format_supported:
        file = os.path.abspath(pathfile)
        logger.debug(f"pathfile: {pathfile}\nfile : {file}")
        if os.path.isfile(file):
            reads = []
            with open(file) as f:
                for read in parse(f, format):
                    logger.info(f"Loading of {read.id}.")
                    reads.append(read)
            return reads
        else:
            raise ValueError("File not found")
    else:
        raise ValueError("Files format not supported")

# TODO Permit the import of folder
