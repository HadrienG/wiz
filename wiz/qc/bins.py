#!/usr/bin/env python
# -*- coding: utf-8 -*-
from wiz.qc.tools import seq_spliter
from wiz.qc import gc, tetra, taxonomy


class Bins:
    def __init__(self, seq, args, gene_pos):
        self.id = seq.id
        self.seq = seq.seq
        self.name = seq.name
        self.gc_percentil = [5, 95]
        self.gc = gc.average_gc(seq_spliter(seq.seq, args.window),
            self.gc_percentil)
        self.gc_global = gc.average_gc_global(self.seq)
        self.gc_bounds = gc.get_bounds(self.gc, self.gc_percentil)
        self.gc_filtered = gc.percentil_filter(self.gc, self.gc_percentil)
        self.tetra = tetra.tetranuc_count(self.seq)
        self.taxonomy = "undefined"
        if self.name in gene_pos.keys():
            self.coding_density = taxonomy.coding_density(gene_pos[self.name],
                len(self.seq))
        else:
            self.coding_density = 0
        

# TODO Validating the class
