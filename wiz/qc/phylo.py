#!/usr/bin/env python
# -*- coding: utf-8

import os
from wiz.qc import tools
import logging


logger = logging.getLogger(__name__)


def phylogeny(filename, args):
    querry = f"{args.output}/prodigal/{filename}.faa"
    output = f"{args.output}/hmmer/hmm_report"
    sequences = {}
    for profil in args.markers:
        tools.hmmscan(querry, profil, output)
        sequences = cleaning_result(output, sequences)
    return sequences


def cleaning_result(outfile, sequences):
    data = []
    with open(outfile, "r") as long_report:
        lines = long_report.readlines()
        for l in lines[3:-10]:
            if "#" not in l:
                data.append(l)
    for line in data:
        l_cut = line.split()
        model, querry, seqEvalue = l_cut[0], l_cut[2], l_cut[4]
        querry = querry.split("_")
        gene = querry[-1]
        querry = querry[:-1]
        querry = ("_".join(querry), gene)
        logger.debug(f""" MATCH!
    Model\t  {model}
    Sequence  {querry[0]}
    Gene\t  {querry[1]}
    Evalue\t  {seqEvalue}""")
        if querry[0] in sequences.keys():
            sequences[querry[0]].append((float(seqEvalue), model, querry[1]))
        else:
            sequences[querry[0]] = [(float(seqEvalue), model, querry[1])]
    return(sequences)