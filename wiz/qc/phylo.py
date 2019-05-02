#!/usr/bin/env python
# -*- coding: utf-8

"""
File that gathers the necessary functions for phylogenics operations.
This function are not used in QC module and must be moved in phylo module.
"""

import logging

from wiz.qc import tools

logger = logging.getLogger(__name__)


def phylogeny(filename, args):
    """
    Compares the output of the prodigal program with a profile
    store via hmmscan
    """
    querry = f"{args.output}/prodigal/{filename}.faa"
    output = f"{args.output}/hmmer/hmm_report"
    sequences = {}
    for profil in args.markers:
        tools.hmmscan(querry, profil, output)
        sequences = cleaning_result(output, sequences)
    return sequences


def cleaning_result(outfile, sequences):
    """
    parse the output file of hmmscan to retrieve the model,
    the sequence and the Evalue
    """
    data = []
    with open(outfile, "r") as long_report:
        lines = long_report.readlines()
        for line in lines[3:-10]:
            if "#" not in line:
                data.append(line)
    for line in data:
        l_cut = line.split()
        model, querry, seq_evalue = l_cut[0], l_cut[2], l_cut[4]
        querry = querry.split("_")
        gene = querry[-1]
        querry = querry[:-1]
        querry = ("_".join(querry), gene)
        logger.debug(f""" MATCH!
    Model\t  {model}
    Sequence  {querry[0]}
    Gene\t  {querry[1]}
    Evalue\t  {seq_evalue}""")
        if querry[0] in sequences.keys():
            sequences[querry[0]].append((float(seq_evalue), model, querry[1]))
        else:
            sequences[querry[0]] = [(float(seq_evalue), model, querry[1])]
    return sequences
