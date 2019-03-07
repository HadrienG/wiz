#!/usr/bin/env python
# -*- coding: utf-8

import logging


logger = logging.getLogger(__name__)


def check_for_duplicates(duplicates_list, auto_filter):
    """Checks based on IDs if duplicates are present"""
    id_list = []
    duplicates_founded = False
    bins_duplicated = []
    for bin in duplicates_list:
        if bin.id not in id_list:
            id_list.append(bin.id)
        else:
            duplicates_founded = True
            bins_duplicated.append(bin.id)
    if duplicates_founded:
        logger.warning("Duplicated bins founded:")
        logger.warn(bins_duplicated)
        if auto_filter:
            logger.info("automatic filtration")
            return automatic_filter(duplicates_list)
        else:
            rep = input("""Do you want :
            \t[D]o nothing
            \t[A]utomatically filter (default)
            =>\t""").lower()
            if rep == "d":
                logger.info("ignored filtration")
                return duplicates_list
            else:
                if not (rep == "a" or rep == ""):
                    logger.warn(
                        f"'{rep}' is a bad value ! The default option is used")
                logger.info("automatic filtration")
                return automatic_filter(duplicates_list)
    else:
        return duplicates_list


def automatic_filter(list_bins):
    """Analyze the list provided and keep the first
    item founded in case of duplicates"""
    id_list = []
    filtered_list = []
    for bin in list_bins:
        if bin.id not in id_list:
            id_list.append(bin.id)
            filtered_list.append(bin)
    return filtered_list


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


def seq_spliter(sequence, window, truncate=True):
    subseqs = []
    if check_window_size(sequence, window):
        average = []
        seq_size = len(sequence)
        for pos in range(0, seq_size, window):
            if pos+window < seq_size:
                subseq = sequence[pos:pos+window]
            else:
                if not truncate:
                    subseq = sequence[pos:]
            subseqs.append(subseq)
    return subseqs
