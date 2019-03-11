#!/usr/bin/env python
# -*- coding: utf-8

from wiz.annotate.tools import prodigal


def coding_density(genes_pos, len_seq):
    coding_region = 0
    overlap_region = 0
    for genes in genes_pos:
        coding_region += (genes[1]-genes[0])
    for genes in genes_pos[:-1]: 
        n_genes = genes_pos[genes_pos.index(genes)+1]
        overlap_region += max(0, min(genes[1],n_genes[1])-max(genes[0],n_genes[0]))
    # print(coding_region,overlap_region)
    # print(((coding_region-overlap_region)/len_seq)*100)
    return ((coding_region-overlap_region)/len_seq)*100


# ok to remove this function
def coding_density2(genes_pos, len_seq):
    # print(genes_pos)
    coding_region = [(genes_pos[0])]
    for genes in genes_pos[1:]:
        r_start, r_stop = coding_region[-1]
        g_start, g_stop = genes
        if g_start > r_stop+1:
            coding_region.append(genes)
        else :
            if g_stop > r_stop :
                coding_region[-1] = (r_start,g_stop)
    tot_coding = 0
    for region in coding_region :
        tot_coding+=(region[1]-region[0])
    # print(tot_coding)
    tot_coding = (tot_coding/len_seq)*100
    # print(tot_coding)
    return tot_coding
# =========================

def taxonomy(parameter_list):
    pass


