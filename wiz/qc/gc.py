#!/usr/bin/env python
# -*- coding:utf-8 -*-


from Bio.SeqUtils import GC
from numpy import percentile

from wiz.qc import tools

import logging
logger = logging.getLogger(__name__)


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


def bounds(sequence, percent=[5, 95]):
    bound_1 = percentile(sequence, percent[0])
    bound_2 = percentile(sequence, percent[1])
    return bound_1, bound_2

