#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
File that gathers the necessary functions for gc content operations.
"""

import logging

from Bio.SeqUtils import GC
from numpy import percentile

from wiz.qc import tools


logger = logging.getLogger(__name__)


def average_gc(sequence, window=-1):
    """
    The default value -1 allows the function to use the size of the
    entire sequence by default. The size of the sequence can be assigned
    directly to the function's call. The transmitted window size is
    automatically checked at the size of the sequence. If the size of
    the window is equal to the size of the sequence the function returns
    a float, if the size of the window is between 0 and the size of the
    sequence, the function returns a list. In other cases an error is thrown.
    """
    average = []
    if window == -1:
        window = len(sequence)
    if tools.check_window_size(sequence, window):
        lenght = len(sequence)
        for pos in range(0, lenght, window):
            if pos+window <= lenght:
                average.append(GC(sequence[pos:pos+window]))
    # ! Warning If the user enters exactly the same size of window
    # ! as the size of the smallest sequence the program may crash
    # TODO fix that
    if window == len(sequence):
        return average[0]
    return average


def bounds(averages, percent=-1):
    """
    Based on a list of averages, the function calculates the X and Y
    percentile bounds defined by the percent = [X, Y] parameter.
    By default X = 5 and Y = 95
    """
    if percent == -1:
        percent = [5, 95]
    bound_1 = percentile(averages, percent[0])
    bound_2 = percentile(averages, percent[1])
    return bound_1, bound_2
