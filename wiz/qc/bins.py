#!/usr/bin/env python
# -*- coding: utf-8 -*-
from tools import seq_spliter

class Bins:
    def __init__(self, seq, window):
        self.id = seq.id
        self.seq = seq.seq
        self.name = seq.name
        self.subseqs = seq_spliter(seq.seq,window)
        self.gc = []
        self.tetra = []
        self.gc_bounds = []
        self.gc_filtered = []
        self.gc_percentil = [5, 95]


# TODO Validating the class
