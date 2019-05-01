#!/usr/bin/env python
# -*- coding: utf-8 -*-

from wiz.qc import tools
from wiz.qc import gc, tetra, taxonomy


class Contig:
    def __init__(self, contig, window):
        self.id = contig.id
        self.name = contig.name
        self.sequence = contig.seq
        # these values are computed during initialisation
        self.length = len(self.sequence)
        self.gc_contig = gc.average_gc(self.sequence)
        self.gc_cutouts = gc.average_gc(self.sequence, window)
        self.gc_bounds = gc.bounds(self.gc_cutouts)
        self.tetranucl = tetra.tetranuc_count(self.sequence)
        # these values are difined after initialisation
        self.coding_density = -25
        self.genes_position = []
        self.taxonomy = [({"species": "undefined", "no rank": "Not found"}, 0)]
        # self.phylogeny_profils = "Not results"


class Bins:
    def __init__(self, path, contigs):
        self.filename = tools.file_name(path)
        self.path = path
        self.contigs = contigs
