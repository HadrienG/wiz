#!/usr/bin/env python
# -*- coding: utf-8
from math import sqrt
from scipy.spatial.distance import cityblock as manhattan
from scipy.spatial.distance import euclidean, chebyshev, minkowski


def tetranuc_count(sequence):
    tetra_dic = {}
    total_count = 0
    buffer = str(sequence[:3])
    for nucl in sequence[3:]:
        buffer += str(nucl)
        total_count +=1
        if buffer in tetra_dic.keys():
            tetra_dic[buffer] += 1
        else:
            tetra_dic[buffer] = 1
        buffer = str(buffer[1:])
    key_list = [key for key in tetra_dic.keys()]
    for key in key_list:
        tetra_dic[key] = (tetra_dic[key]/total_count)
    return tetra_dic


# * It's ok here
# TODO Make dataframe with the values returned
def tetra_euclidian_distance(seq1={}, seq2={}):
    # https://en.wikipedia.org/wiki/Euclidean_distance
    dimensions = [key for key in seq1.keys()]
    for key in seq2.keys():
        if key not in dimensions:
            dimensions.append(key)
    total = 0
    for dim in dimensions:
        val_seq1, val_seq2 = 0, 0
        if dim in seq1.keys():
            val_seq1 = seq1[dim]
        if dim in seq2.keys():
            val_seq2 = seq2[dim]
        total += (val_seq1-val_seq2)**2
    return sqrt(total)


def tetra_manhattan_distance(seq1={}, seq2={}):
    # https://fr.wikipedia.org/wiki/Distance_de_Manhattan
    dimensions = [key for key in seq1.keys()]
    for key in seq2.keys():
        if key not in dimensions:
            dimensions.append(key)
    total = 0
    for dim in dimensions:
        val_seq1, val_seq2 = 0, 0
        if dim in seq1.keys():
            val_seq1 = seq1[dim]
        if dim in seq2.keys():
            val_seq2 = seq2[dim]
        total += abs(val_seq2-val_seq1)
    return total


def tetra_scipy_manhattan(seq1={}, seq2={}):
    dimensions = [key for key in seq1.keys()]
    for key in seq2.keys():
        if key not in dimensions:
            dimensions.append(key)
    val_seq1, val_seq2 = [], []
    for dim in dimensions:
        val1, val2 = 0, 0
        if dim in seq1.keys():
            val1 = seq1[dim]
        if dim in seq2.keys():
            val2 = seq2[dim]
        val_seq1.append(val1)
        val_seq2.append(val2)
    return manhattan(val_seq1, val_seq2)


def tetra_scipy_euclidian(seq1={}, seq2={}):
    dimensions = [key for key in seq1.keys()]
    for key in seq2.keys():
        if key not in dimensions:
            dimensions.append(key)
    val_seq1, val_seq2 = [], []
    for dim in dimensions:
        val1, val2 = 0, 0
        if dim in seq1.keys():
            val1 = seq1[dim]
        if dim in seq2.keys():
            val2 = seq2[dim]
        val_seq1.append(val1)
        val_seq2.append(val2)
    return euclidean(val_seq1, val_seq2)


def tetra_scipy_chebyshev(seq1={}, seq2={}):
    dimensions = [key for key in seq1.keys()]
    for key in seq2.keys():
        if key not in dimensions:
            dimensions.append(key)
    val_seq1, val_seq2 = [], []
    for dim in dimensions:
        val1, val2 = 0, 0
        if dim in seq1.keys():
            val1 = seq1[dim]
        if dim in seq2.keys():
            val2 = seq2[dim]
        val_seq1.append(val1)
        val_seq2.append(val2)
    return chebyshev(val_seq1, val_seq2)


def tetra_scipy_minkowski(seq1={}, seq2={}):
    dimensions = [key for key in seq1.keys()]
    for key in seq2.keys():
        if key not in dimensions:
            dimensions.append(key)
    val_seq1, val_seq2 = [], []
    for dim in dimensions:
        val1, val2 = 0, 0
        if dim in seq1.keys():
            val1 = seq1[dim]
        if dim in seq2.keys():
            val2 = seq2[dim]
        val_seq1.append(val1)
        val_seq2.append(val2)
    return minkowski(val_seq1, val_seq2)


def tetra_distance(bins):
    distances = {}
    for binI in bins[:-1]:
        for binJ in bins[bins.index(binI)+1:]:
            distances[(binI.id, binJ.id)] = {"manhattan":tetra_manhattan_distance(binI.tetra, binJ.tetra)} #,"manhattan":tetra_scipy_manhattan(binI.tetra, binJ.tetra),"chebyshev":tetra_scipy_chebyshev(binI.tetra, binJ.tetra)}
    return distances
