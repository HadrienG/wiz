#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from wiz.qc import tetra


class Contig:
    def __init__(self, contig, window):
        self.id = contig.id
        self.name = contig.name
        self.sequence = contig.seq

path = [
    "wiz/test/Escherichia_coli.fasta",
    "wiz/test/saccharomyces_cerevisiae_CH_II.fasta"
    ]
contigs = []

for p in path:
    f = open(p, "r")
    with f:
        fasta_file = SeqIO.parse(f, "fasta")
        for record in fasta_file:
            contigs.append(Contig(record, 10))

value = {('NC_000913.3', '>tpg|BK006936.2|'): 0.468044738948156}
assert value == tetra.distance_calculation(contigs)
