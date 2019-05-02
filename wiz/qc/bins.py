#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File that contains the objects used in the qc module.
"""

from wiz.qc import gc
from wiz.qc import tetra
from wiz.qc import tools


class Contig:
    """
        The contig object represents a contig with its characteristics
        (name, id, sequence), metrics made at its creation (length,% gc,
        tetranucleic composition, ..) and initialized metrics with defaults
        to be attributed later because more consistent in calculations.
        """
    def __init__(self, contig, window):
        self.uid = contig.id
        self.name = contig.name
        self.sequence = contig.seq
        # these values are computed during initialisation
        self.length = len(self.sequence)
        self.gc_contig = gc.average_gc(self.sequence, self.length)
        self.gc_cutouts = gc.average_gc(self.sequence, window)
        self.gc_bounds = gc.bounds(self.gc_cutouts)
        self.tetranucl = tetra.tetranuc_count(self.sequence)
        # these variables are defined after initialization
        self.coding_density = -25
        self.genes_position = []
        self.taxonomy = [({"species": "undefined", "no rank": "Not found"}, 0)]
        # self.phylogeny_profils = "Not results"


class Bins:
    """
    The bin object represents a set of contigs from the binning,
    the name of the original file and its path.
    """
    def __init__(self, path, contigs):
        self.filename = tools.file_name(path)
        self.path = path
        # The contigs variable is a contig object list
        self.contigs = contigs
