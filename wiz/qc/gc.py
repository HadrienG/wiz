#!/usr/bin/env python
# -*- coding:utf-8 -*-


from Bio.SeqUtils import GC
from numpy import percentile

from wiz.qc import tools

import logging
logger = logging.getLogger(__name__)


# def average_gc(subseqs, size=len(seq) truncate=False):  # I think it's OK
#     """
#     Return a list with for each window of size window_size the average of GC.
#     The last one windows can be more small that the others windows due of the
#     sequence length. The result of the window of length inferior at 
#     window_size can be removed with the parameter truncate=True
#     """
#     average = []
#     for subseq in subseqs:
#         average.append(GC(subseq))
#     return average


# def average_gc_global(sequence):
#     return(GC(sequence))

def average_gc(sequence, window=-1):
    average = []
    if window == -1:
        window = len(sequence)
    if tools.check_window_size(sequence, window):
        lenght = len(sequence)
        for pos in range(0, lenght, window):
            if pos+window <= lenght:
                average.append(GC(sequence[pos:pos+window]))
    if len(average) == 1:
        return average[0]
    else:
        return average


# def percentil_filter(average, percent=[5, 95]):  # WIP need testing
#     P5, P95 = get_bounds(average, percent)
#     # logger.debug(f"\t| P5 : {round(P5,3):>6} |\tP95 : {round(P95,3):>6}|")
#     filtered_average = []
#     for value in average:
#         if value < P95 and value > P5:
#             filtered_average.append(value)
#     return filtered_average


def bounds(sequence, percent=[5, 95]):
    bound_1 = percentile(sequence, percent[0])
    bound_2 = percentile(sequence, percent[1])
    return bound_1, bound_2


# TODO Try to plot on graph windows outside the 95e and the 5e percentil
# * Data filtered before that
