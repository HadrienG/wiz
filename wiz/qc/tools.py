#!/usr/bin/env python
# -*- coding: utf-8

from Bio.SeqIO import parse
import os
import logging


logger = logging.getLogger(__name__)


def genome_parser(pathfiles, format="fasta"):
    logger.debug("[START] genome_parser")
    reads = []
    genpath = os.path.abspath(pathfiles)
    if os.path.isfile(genpath):
        logger.debug("genome_parser : file detected")
        reads.append(file_reader(genpath, format)[0])
    elif os.path.isdir(genpath):
        logger.debug("genome_parser : folder detected")
        reads.append(folder_reader(genpath, format)[0])
    else:
        raise ValueError("File not found")
    logger.debug("[STOP] genome_parser")
    return reads


def file_reader(file_path, format):
    logger.debug(f"Opening file : {file_path}")
    reads = []
    with open(file_path) as f:
        for read in parse(f, format):
            logger.info(f"Loading of {read.id}.")
            reads.append(read)
        return reads


def folder_reader(folder_path, format):
    logger.debug("[START] folder_reader")
    file_list = os.listdir(folder_path)
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
                    reads.append(file_reader(file_path, format)[0])
            else:
                file_path = os.path.join(folder_path, file_list[0])
                reads.append(file_reader(file_path, format)[0])
            logger.debug("[STOP] folder_reader")
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
