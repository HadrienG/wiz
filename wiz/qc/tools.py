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
            if pos+window < seq_size:
                subseq = sequence[pos:pos+window]
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
    with open(f"{os.path.abspath(path)}/{file}.gff",'r') as f:
        for i in f.readlines():
            if not "#" in i:
                cut=i.split("\t")
                if cut[0] in sequence.keys():
                    sequence[cut[0]].append((int(cut[3]),int(cut[4])))
                else:
                    sequence[cut[0]] = [(int(cut[3]),int(cut[4]))]
    return sequence


def finch(genome_id, path_db, output_dir):
    finch = path.software_exists("finch")
    logger.info(" Running Finch")
    logger.info(f" sketching {genome_id}")
    subprocess.run([finch, "sketch",f"{output_dir}/finch/{genome_id}.fna","-o", f"{output_dir}/finch/{genome_id}.sk"])
    logger.info(f" Compare {genome_id} to DB")
    subprocess.run([finch,"dist","-o", f"{output_dir}/finch/{genome_id}.finchout", path_db, f"{output_dir}/finch/{genome_id}.sk"])#, capture_output=True)
    #with open(f"{output_dir}/finch/{genome_id}.finchout","ab") as f:
    #    f.write(output.stdout)

