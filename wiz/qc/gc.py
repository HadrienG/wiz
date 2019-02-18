#!/usr/bin/env python
# -*- coding:utf-8 -*-


from Bio.SeqUtils import GC
from numpy import percentile


def check_window_size(sequence, window_size):       # I think it's OK
    if len(sequence) == 0:
        error = "The sequence is void."
        raise ValueError(error)
    if window_size < 0:
        error = "The size of the window is negative."
        raise ValueError(error)
    if window_size == 0:
        error = "The size of the window is null."
        raise ValueError(error)
    if window_size > len(sequence):
        error = "The size of the window is superior of the sequence length."
        raise ValueError(error)
    return True


def average_gc(sequence, window_size=5000, truncate=False):  # I think it's OK
    """
    Return a list with for each window of size window_size the average of GC.
    The last one windows can be more small that the others windows due of the
    sequence length. The result of the window of length inferior at window_size
    can be removed with the parameter truncate=True
    """
    if check_window_size(sequence, window_size):
        average = []
        seq_size = len(sequence)
        for pos in range(0, seq_size, window_size):
            if pos+window_size < seq_size:
                subseq = sequence[pos:pos+window_size]
            else:
                if not truncate:
                    subseq = sequence[pos:]
                else:
                    break
            average.append(GC(subseq))
        return average


def percentil_filter(average, percent=[5, 95]):  # WIP need testing
    P5, P95 = percentile(average, percent[0]), percentile(average, percent[1])
    for value in average:
        if value > P95 or value < P5:
            average.pop(average.index(value))
    return average

# TODO Try to plot on graph windows outside the 95e and the 5e percentil
# * Data filtered before that
