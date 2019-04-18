#!/usr/bin/env python
# -*- coding: utf-8
<<<<<<< HEAD
import os
from wiz.qc.tools import hmmscan
import logging


logger = logging.getLogger(__name__)


def phylogeny(genome, args):
    file = os.path.basename(genome)
    if "." in file:
        file = file.split(".")[0]
    results = targeting(
        f"{args.output}/prodigal/{file}.faa",
        f"{args.hmmdb}targeting/",
        f'{args.output}/hmmer/')
    return results


def targeting(querry, db_folder, outdir):
    dbs = os.listdir(db_folder)
    hmm_profils = []
    first_pass = True
    sequences = {}
    for db in dbs:
        if not".hmm." in db:
            hmm_profils.append(db)
    for profil in hmm_profils:
        logger.debug(f"hmmscan with : {db_folder}{profil}")
        hmmscan(querry, f"{db_folder}{profil}", f"{outdir}hmm_report")
        sequences = cleaning_result(outdir, first_pass, sequences)
        first_pass = False
    return(sequences)


def identification():
    pass


def cleaning_result(outdir, header, sequences):
    to_write = []
    header_text = []
    with open(f"{outdir}hmm_report", "r") as long_report:
        lines = long_report.readlines()
        header_text = lines[:3]
        for l in lines[:-10]:
            if "#" not in l:
                to_write.append(l)
    if header:
        with open(f"{outdir}hmm_result", "w") as short_report:
            for l in header_text:
                short_report.writelines(l)
            if len(to_write) > 0:
                for l in to_write:
                    short_report.writelines(l)
    with open(f"{outdir}hmm_result", "a") as short_report:
        if len(to_write) > 0:
            for l in to_write:
                short_report.writelines(l)
    for l in to_write:
        l_cut = l.split()
        logger.debug(l_cut)
        model, querry, seqEvalue = l_cut[0], l_cut[2], l_cut[4]
        querry = querry.split("_")
        gene = querry[-1]
        querry = querry[:-1]
        querry = ("_".join(querry), gene)
        logger.debug(f"{model}/seq {querry[0]}g{querry[1]}: seq Evalue {seqEvalue}")
        if querry[0] in sequences.keys():
            sequences[querry[0]].append((float(seqEvalue), model, querry[1]))
        else:
            sequences[querry[0]] = [(float(seqEvalue), model, querry[1])]
    return(sequences)
