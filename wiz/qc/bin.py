#!/usr/bin/env python
# -*- coding: utf-8 -*-
from wiz.qc.tools import seq_spliter
from wiz.qc import gc, tetra


class Bin:
    def __init__(self, seq, window):
        self.id = seq.id
        self.seq = seq.seq
        self.name = seq.name
        self.subseqs = seq_spliter(seq.seq, window)
        self.gc_percentil = [5, 95]
        self.gc = gc.average_gc(self.subseqs, self.gc_percentil)
        self.gc_bounds = gc.get_bounds(self.gc, self.gc_percentil)
        self.gc_filtered = gc.percentil_filter(self.gc, self.gc_percentil)
        self.tetra = tetra.tetranuc_count(self.seq)


# TODO Validating the class
