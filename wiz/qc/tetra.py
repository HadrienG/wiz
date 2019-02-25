#!/usr/bin/env python
# -*- coding: utf-8
from math import sqrt
from scipy.spatial.distance import cityblock as manhattan
from scipy.spatial.distance import euclidean, chebyshev, minkowski
import threading
import queue
from os import cpu_count
from time import sleep

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


# * It's ok here
# TODO Make dataframe with the values returned
# def tetra_euclidian_distance(seq1={}, seq2={}):
#     # https://en.wikipedia.org/wiki/Euclidean_distance
#     dimensions = [key for key in seq1.keys()]
#     for key in seq2.keys():
#         if key not in dimensions:
#             dimensions.append(key)
#     total = 0
#     for dim in dimensions:
#         val_seq1, val_seq2 = 0, 0
#         if dim in seq1.keys():
#             val_seq1 = seq1[dim]
#         if dim in seq2.keys():
#             val_seq2 = seq2[dim]
#         total += (val_seq1-val_seq2)**2
#     return sqrt(total)


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


# def tetra_scipy_manhattan(seq1={}, seq2={}):
#     dimensions = [key for key in seq1.keys()]
#     for key in seq2.keys():
#         if key not in dimensions:
#             dimensions.append(key)
#     val_seq1, val_seq2 = [], []
#     for dim in dimensions:
#         if dim in seq1.keys():
#             val_seq1.append(seq1[dim])
#         else:
#             val_seq1.append(0)
#         if dim in seq2.keys():
#             val_seq2.append(seq2[dim])
#         else:
#             val_seq2.append(0)
#     return manhattan(val_seq1, val_seq2)


# def tetra_scipy_euclidian(seq1={}, seq2={}):
#     dimensions = [key for key in seq1.keys()]
#     for key in seq2.keys():
#         if key not in dimensions:
#             dimensions.append(key)
#     val_seq1, val_seq2 = [], []
#     for dim in dimensions:
#         val1, val2 = 0, 0
#         if dim in seq1.keys():
#             val1 = seq1[dim]
#         if dim in seq2.keys():
#             val2 = seq2[dim]
#         val_seq1.append(val1)
#         val_seq2.append(val2)
#     return euclidean(val_seq1, val_seq2)


# def tetra_scipy_chebyshev(seq1={}, seq2={}):
#     dimensions = [key for key in seq1.keys()]
#     for key in seq2.keys():
#         if key not in dimensions:
#             dimensions.append(key)
#     val_seq1, val_seq2 = [], []
#     for dim in dimensions:
#         val1, val2 = 0, 0
#         if dim in seq1.keys():
#             val1 = seq1[dim]
#         if dim in seq2.keys():
#             val2 = seq2[dim]
#         val_seq1.append(val1)
#         val_seq2.append(val2)
#     return chebyshev(val_seq1, val_seq2)


# def tetra_scipy_minkowski(seq1={}, seq2={}):
#     dimensions = [key for key in seq1.keys()]
#     for key in seq2.keys():
#         if key not in dimensions:
#             dimensions.append(key)
#     val_seq1, val_seq2 = [], []
#     for dim in dimensions:
#         val1, val2 = 0, 0
#         if dim in seq1.keys():
#             val1 = seq1[dim]
#         if dim in seq2.keys():
#             val2 = seq2[dim]
#         val_seq1.append(val1)
#         val_seq2.append(val2)
#     return minkowski(val_seq1, val_seq2)


# def tetra_distance(bins):
    
#     distances = []
#     # task = ((len(bins)**2)-len(bins))/2
#     # count = 0
#     for binI in bins[:-1]:
#         for binJ in bins[bins.index(binI)+1:]:
#             tasks.append(((binI.id, binJ.id), (binI.tetra, binJ.tetra)))
#     print("para")
#     for task in tasks[:1]:
#         print(task)
#     jobs = [distances_calc(tasks, distances) for i in range(0, cpu_count()-1)]
#     for job in jobs:
#         job.start()
#     for job in jobs:
#         job.join()
#     #Parallel(n_jobs=2)((task[0], tetra_manhattan_distance(task[1][0], task[1][1])) for task in tasks)
#     # for task in tasks:
#         # distances[task[0]] = tetra_manhattan_distance(task[1][0], task[1][1])
#     return distances

exitFlag = 0
queueLock = threading.Lock()
workQueue = queue.Queue(0)
tetra_dist = {}


class Tdistances_calc(threading.Thread):
    def __init__(self, ID, queue):
        threading.Thread.__init__(self)
        self.id = ID
        self.q = queue

    def run(self):
        logger.info(f"Starting Thread #{self.id}")
        process_data(self.id, self.q)
        logger.info(f"Exiting Thread #{self.id}")


def process_data(thread_id, queue):
    dead_count = 0
    while not exitFlag:
        queueLock.acquire()
        if not workQueue.empty():
            id, task = queue.get()
            logger.debug(f"Thread {thread_id} task {id}")
            queueLock.release()
            result = tetra_manhattan_distance(task[0], task[1])
            queueLock.acquire()
            tetra_dist[id] = result
            queueLock.release()
        else:
            logger.debug(f"Thread {thread_id} waiting")
            queueLock.release()
            dead_count += 1
            if dead_count < 5:
                sleep(1)
            else:
                break


def distance_calculation(bins):
    threads = []
    threadID = 1

    tasks = []
    for binI in bins[:-1]:
        for binJ in bins[bins.index(binI)+1:]:
            tasks.append(((binI.id, binJ.id), (binI.tetra, binJ.tetra)))
    tasks_nb = len(tasks)

    for i in range(0, cpu_count()*2):
        thread = Tdistances_calc(threadID, workQueue)
        thread.start()
        threads.append(thread)
        threadID += 1
    
    queueLock.acquire()
    for task in tasks:
        workQueue.put(task)
    queueLock.release()

    while True:
        logger.info(f"completed tasks :\t{tasks_nb-workQueue.qsize()} \t/{tasks_nb}")
        if workQueue.qsize() == 0:
            break
        else:
            sleep(1)
            pass

    exitFlag = 1

    for t in threads:
        t.join()

    logger.info("Exiting tetranucleic distance calculation module")
    return tetra_dist