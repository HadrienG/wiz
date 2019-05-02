#!/usr/bin/env python
# -*- coding: utf-8

"""
File that gathers the necessary functions for tetranucleic operations.
"""

from os import cpu_count
from multiprocessing import Pool


def tetranuc_count(sequence):
    """
    Function who's calculating the proportion of each
    tetranucleotid in a sequence.
    """
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


def tetra_manhattan_distance(seq1, seq2):
    """
    A function that calculates the Manhattan distance between two
    dictionaries representative of the proportion of tetranucleotides
    in each sequence.
    https://fr.wikipedia.org/wiki/Distance_de_Manhattan
    """
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


def distance_calculation(contigs, nb_process=(cpu_count()-1)):
    """
    Function that prepares the tasks for the pool calculates
    and launches this one
    """
    tasks = []
    for contig_i in contigs[:-1]:
        for contig_j in contigs[contigs.index(contig_i)+1:]:
            tasks.append(
                (
                    (
                        contig_i.uid,
                        contig_j.uid
                    ),
                    (
                        tetranuc_count(contig_i.sequence),
                        tetranuc_count(contig_j.sequence)
                    )
                ))
    results = []
    with Pool(processes=nb_process) as pool:
        results = pool.map(compute_dist, tasks)
    dict_results = {}
    for result in results:
        contig_id, value = result
        dict_results[contig_id] = value
    return dict_results


def compute_dist(task):
    """
    Function used in the parallelization of Manhattan distance
    calculation tasks between two sequences.
    """
    id_task, values = task
    result = tetra_manhattan_distance(values[0], values[1])
    return (id_task, result)
