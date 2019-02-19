#!/usr/bin/env python
# -*- coding: utf-8

from Bio.SeqIO import parse
import os
import logging


logger = logging.getLogger(__name__)


def genome_parser(pathfiles, format="fasta"):
    """Checks the supplied paths and calls reads retrieval or folder open functions"""
    logger.debug("[START] genome_parser")
    reads = []
    genpath = os.path.abspath(pathfiles)
    if os.path.isfile(genpath):
        logger.debug("genome_parser : file detected")
        reads += file_reader(genpath, format)
    elif os.path.isdir(genpath):
        logger.debug("genome_parser : folder detected")
        reads += folder_reader(genpath, format)
    else:
        raise ValueError("File not found")
    logger.debug("[STOP] genome_parser")
    return reads


def file_reader(file_path, format):
    """Extract the reads from the provided files"""
    logger.debug(f"Opening file : {file_path}")
    reads = []
    with open(file_path) as f:
        for read in parse(f, format):
            logger.info(f"Loading of {read.id}.")
            reads.append(read)
        return reads


def folder_reader(folder_path, format):
    """List the files in the folders and propose their extractions"""
    logger.debug("[START] folder_reader")
    file_list = os.listdir(folder_path)
    file_list_temps = []
    for file in file_list:
        path_f = os.path.join(folder_path, file)
        if os.path.isfile(path_f):
            file_list_temps.append(file)
    file_list = list(file_list_temps)
    del(file_list_temps)
    folder_name = os.path.split(folder_path)[1]
    logger.warn("The provided path points to a folder.")
    if len(file_list) != 0:
        if len(file_list) == 1:
            end_string = "only one file"
        else:
            end_string = f"these {len(file_list)} files"
        rep = input(f"Do you want to open the folder '{folder_name}' with {end_string} ? [y/N]")
        if rep == "y":
            reads = []
            logger.debug(f"lengt of file list : {len(file_list)}")
            if len(file_list) > 1:
                for file in file_list:
                    file_path = os.path.join(folder_path, file)
                    reads += file_reader(file_path, format)
            else:
                file_path = os.path.join(folder_path, file_list[0])
                reads += file_reader(file_path, format)
            logger.debug("[STOP] folder_reader")
            logger.debug(f"Folder_reader : reads => {reads}")
            return reads
        else:
            logger.warning("Interrupted operations")
            logger.debug("[STOP] folder_reader")
            exit(0)
    else:
        logger.warning(f"The folder '{folder_name}' is empty")
        logger.debug("[STOP] folder_reader")
        exit(0)
# * OK >> TODO Permit the import of folder


def check_for_duplicates(duplicates_list):
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
        rep = input("""Do you want :
        \t[D]o nothing
        \t[A]utomatically filter (default)
        \t[M]anually filter
        =>\t""").lower()
        if rep == "d":
            logger.info("ignored filtration")
            return duplicates_list
        elif rep == "m":
            return manually_filter(duplicates_list, bins_duplicated)
        else:
            if not (rep == "a" or rep == ""):
                logger.warn(f"'{rep}' is a bad value ! The default option is used")
            logger.info("automatic filtration")
            return automatic_filter(duplicates_list)
    else:
        return duplicates_list


def automatic_filter(list_bins):
    """Analyze the list provided and keep the first item found in case of duplicates"""
    id_list = []
    filtered_list = []
    for bin in list_bins:
        if bin.id not in id_list:
            id_list.append(bin.id)
            filtered_list.append(bin)
    return filtered_list


def manually_filter(list_bins, bins_duplicated):
    logger.warning("Manual filtration")
    logger.info(f"{len(list_bins)} bin involved")
    list_bins.reverse()
    pos_bin = [bin.id for bin in list_bins]
    for dupl in bins_duplicated:
        rep = input(f"{dupl} : keep ? [y/N]").lower()
        if not rep == "y":
            list_bins.pop(pos_bin.index(dupl))
            pos_bin.pop(pos_bin.index(dupl))
            logger.info(f"{dupl} removed")
    list_bins.reverse()
    return list_bins
