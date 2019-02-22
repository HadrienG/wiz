#!/usr/bin/env python
# -*- coding: utf-8 -*-


class Bin:
    def __init__(self, seq):
        self.id = seq.id
        self.seq = seq.seq
        self.name = seq.name
        self.gc = []
        self.tetra = []
        self.gc_bounds = []
        self.gc_filtered = []
        self.gc_percentil = [5, 95]


# TODO Validating the class
