#!/usr/bin/env python
# -*- coding: utf-8 -*-
from wiz.qc.tools import seq_spliter
from wiz.qc import gc, tetra, taxonomy


class Bins:
    def __init__(self, seq, args, gene_pos, filename, pos_in_file, phylo):
        self.origin_file = filename
        self.pos_in_file = pos_in_file
        self.id = seq.id
        self.seq = seq.seq
        self.name = seq.name
        self.gc_percentil = [5, 95]
        self.gc = gc.average_gc(
            seq_spliter(seq.seq, args.window),
            self.gc_percentil)
        self.gc_global = gc.average_gc_global(self.seq)
        self.gc_bounds = gc.get_bounds(self.gc, self.gc_percentil)
        self.gc_filtered = gc.percentil_filter(self.gc, self.gc_percentil)
        self.tetra = tetra.tetranuc_count(self.seq)
        self.taxonomy = [({"species": "undefined", "no rank": "Not found"}, 0)]
        if self.name in gene_pos.keys():
            self.coding_density = taxonomy.coding_density(
                gene_pos[self.name],
                len(self.seq))
        elif f"Prodigal_Seq_{str(self.pos_in_file)}" in gene_pos.keys():
            self.coding_density = taxonomy.coding_density(
                gene_pos[f"Prodigal_Seq_{str(self.pos_in_file)}"],
                len(self.seq))
        else:
            self.coding_density = -25
        if self.name in phylo.keys():
            self.phylogeny_profils = phylo[self.name]
        elif f"Prodigal_Seq_{self.pos_in_file}" in phylo.keys():
            self.phylogeny_profils = phylo[f"Prodigal_Seq_{self.pos_in_file}"]
        else:
            self.phylogeny_profils = "No results"


# TODO Validating the class
