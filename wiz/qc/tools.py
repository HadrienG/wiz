#!/usr/bin/env python
# -*- coding: utf-8

"""
File that gathers the necessary functions and tools for wiz qc.
"""

import os
import logging
import subprocess

from wiz.misc import path

logger = logging.getLogger(__name__)


def create_dirs(output, force_overwrite):
    """
    Function that creates the folders needed by the program.
    """
    logger.debug(" Creating dirs")
    if not os.path.exists(output) or force_overwrite:
        if force_overwrite:
            logger.info(" You used the --force young padawan!")
            if os.path.exists(output):
                logger.warning(" The output folder already exists!")
                logger.warning(" The data in it may be overwritten!")
        path.create_dir(output, force=True)
        path.create_dir(output+"/prodigal", force=True)
        path.create_dir(output+"/hmmer", force=True)
        path.create_dir(output+"/finch", force=True)
        logger.info(" The path folder as been created.")
        logger.info(f" All output files will be stored in {output}")
    else:
        logger.warning(f" The path {output} already exists !")
        logger.info(" use --force to overwrite the data")
        exit(2)  # error code for "command line error in unix like system"


def check_window_size(sequence, window_size):
    """
    check the value of the window size and the length of the sequence
    the window_size must be 0 < size <= lenght(sequence)
    the lenght must be 0 < lenght(sequence)
    """
    if not sequence:  # len(sequence) == 0:
        error = "The sequence is void."
        raise ValueError(error)
    if window_size < 0:
        error = "The size of the window is negative."
        raise ValueError(error)
    if window_size == 0:
        error = "The size of the window is null."
        raise ValueError(error)
    if window_size > len(sequence):
        error = f"The size of the window ({window_size}) is "
        error += f"superior of the sequence length ({len(sequence)})."
        raise ValueError(error)
    return True


def file_name(file):
    """
    function who extract the name in a path
    """
    file = os.path.basename(file)
    # if filename = file.ext the part ext is removed
    # ! problem if filename is file.10 because the number is removed
    if '.' in file:
        file = file.split(".")[:-1]
        file = ".".join(file)
    return file


def genes(path_file, file):
    """
    Function who's extracting the genes positions inside a
    prodigal output file.
    """
    sequence = {}
    with open(f"{os.path.abspath(path_file)}/{file}.gff", 'r') as read_file:
        for i in read_file.readlines():
            # if the line is not a Prodigal comment
            if "#" not in i:
                cut = i.split("\t")
                # cut[0] = contig id
                # cut[3] = gene start position
                # cut[4] = gene stop position
                if cut[0] in sequence.keys():
                    # if dictionnary "sequence" contain
                    # already the key contig_id
                    sequence[cut[0]].append((int(cut[3]), int(cut[4])))
                else:
                    # if dictionnary "sequence" don't contain the key contig_id
                    sequence[cut[0]] = [(int(cut[3]), int(cut[4]))]
    return sequence


def finch_sketch(filename, output):
    """
    prepare a sketch of the "filename" to being compared with a sketchsbank
    """
    finch = path.software_exists("finch")
    logger.info(f" Finch sketching contigs in {filename}")
    s_args = [
        finch,
        "sketch",
        "-N",
        "-o",
        f"{output}/finch/{filename}.sk",
        f"{output}/finch/{filename}.fna"]
    subprocess.run(s_args)
    logger.debug(" Finch sketching end")


def finch_dist(filename, dbs, outdir):
    """
    search distance between sketchsbank and a sketch of file "filename"
    """
    logger.info(f" Finch compare {filename} to DB")
    finch = path.software_exists("finch")
    s_args = [
        finch,
        "dist",
        "--max-dist",
        "0.2",
        "-o",
        f"{outdir}/finch/{filename}"]
    s_args += dbs
    s_args += [f"{outdir}/finch/{filename}.sk"]
    subprocess.run(s_args)
    logger.debug(" Finch conparating end")


def hmmscan(filename, profil, outdir):
    """
    Make a hmmscan research
    """
    txt = [file_name(filename), file_name(profil)]
    logger.info(f" Hmmer compare {txt[0]} and {txt[1]}")
    # compare a marker profil and proteins in file "filename"
    hmmer = path.software_exists("hmmscan")
    s_args = [
        hmmer,
        "--noali",
        # "--max",
        "--tblout",
        outdir,
        profil,
        filename
    ]
    fnull = open(os.devnull, "w")  # opening /dev/null like a file
    subprocess.run(s_args, stdout=fnull)  # writing stdout in /dev/null
    logger.debug(" HMMScan end")
