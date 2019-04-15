#!/usr/bin/env python
# -*- coding: utf-8

import logging
import os
from wiz.misc import path
import subprocess

logger = logging.getLogger(__name__)


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
            # subseq = ""
            if pos+window <= seq_size:
                subseq = sequence[pos:pos+window]
                subseqs.append(subseq)
            else:
                if not truncate:
                    subseq = sequence[pos:]
                    subseqs.append(subseq)
    return subseqs


def get_file_name(file):
    file = os.path.basename(file)
    if '.' in file:
        file = file.split(".")[0]
    return file


def get_gene_seq(path, file):
    sequence = {}
    with open(f"{os.path.abspath(path)}/{file}.gff", 'r') as f:
        for i in f.readlines():
            if "#" not in i:
                cut = i.split("\t")
                if cut[0] in sequence.keys():
                    sequence[cut[0]].append((int(cut[3]), int(cut[4])))
                else:
                    sequence[cut[0]] = [(int(cut[3]), int(cut[4]))]
    return sequence


def finch(genome_id, path_db, output_dir):
    finch = path.software_exists("finch")
    logger.debug(" Running Finch")
    logger.debug(f" sketching {genome_id}")
    subprocess.run([
        finch,
        "sketch",
        f"{output_dir}/finch/{genome_id}.fna",
        "-o",
        f"{output_dir}/finch/{genome_id}.sk"])
    logger.debug(f" Compare {genome_id} to DB")
    s_args = [
        finch,
        "dist",
        "--max-dist",
        "0.2",
        "-o",
        f"{output_dir}/finch/{genome_id}"]
    s_args += path_db
    s_args += ["--queries", f"{output_dir}/finch/{genome_id}.sk"]
    subprocess.run(s_args)


def finch_sketch(filename, output):
    finch = path.software_exists("finch")
    logger.debug(f" Finch sketching contigs in {filename}")
    subprocess.run([
        finch,
        "sketch",
        "-N",
        "-o",
        f"{output}/finch/{filename}.sk",
        f"{output}/finch/{filename}.fna"])
    logger.debug(" Finch sketching end")


def finch_dist(filename, dbs, output_dir):
    logger.debug(f" Finch Compare {filename} to DB")
    finch = path.software_exists("finch")
    s_args = [
        finch,
        "dist",
        "--max-dist",
        "0.2",
        "-o",
        f"{output_dir}/finch/{filename}"]
    s_args += dbs
    s_args += [f"{output_dir}/finch/{filename}.sk"]
    subprocess.run(s_args)
    logger.debug(" Finch conparating end")
