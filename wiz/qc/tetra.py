#!/usr/bin/env python
# -*- coding: utf-8
from os import cpu_count
from multiprocessing import Pool

import logging
logger = logging.getLogger(__name__)


def tetranuc_count(sequence):
    tetra_dic = {}
    total_count = len(sequence)-3
    buffer = str(sequence[:3])
    for nucl in sequence[3:]:
        buffer += str(nucl)
        if buffer in tetra_dic.keys():
            tetra_dic[buffer] += 1
        else:
            tetra_dic[buffer] = 1
        buffer = str(buffer[1:])
    key_list = [key for key in tetra_dic.keys()]
    for key in key_list:
        tetra_dic[key] = (tetra_dic[key]/total_count)
    return tetra_dic


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


def distance_calculation(bins, nb_process=(cpu_count()-1)):
    tasks = []
    for binI in bins[:-1]:
        for binJ in bins[bins.index(binI)+1:]:
            tasks.append(((binI.id, binJ.id), (binI.tetra, binJ.tetra)))
    results = []
    with Pool(processes=nb_process) as pool:
        results = pool.map(SPcalculation, tasks)
    dict_results = {}
    for res in results:
        id, value = res
        dict_results[id] = value
    return dict_results


def SPcalculation(task):
    id, values = task
    result = tetra_manhattan_distance(values[0], values[1])
    return (id, result)
