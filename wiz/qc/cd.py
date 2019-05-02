#!/usr/bin/env python
# -*- coding: utf-8

"""
File that gathers the necessary functions for coding density operations.
"""

import logging

logger = logging.getLogger(__name__)


def coding_density(genes_pos, len_seq):
    """
    Function that calculates the coding density of a sequence based on
    its length and the position of its genes.
    """
    coding_region = 0
    overlap_region = 0
    for genes in genes_pos:
        coding_region += (genes[1]-genes[0])
    for genes in genes_pos[:-1]:
        n_genes = genes_pos[genes_pos.index(genes)+1]
        overlap_region += max(
            0,
            min(genes[1], n_genes[1])-max(genes[0], n_genes[0]))
    coding = ((coding_region-overlap_region)/len_seq)*100
    logger.debug(f" Coding density : {coding}")
    return ((coding_region-overlap_region)/len_seq)*100
