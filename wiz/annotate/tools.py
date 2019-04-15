#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import subprocess

from wiz.misc import path


def prodigal(genome, output_prefix="bin", trans_table=11, output_dir=None):
    """wrapper function around prodigal
    """
    logger = logging.getLogger(__name__)

    prodigal = path.software_exists("prodigal")
    if output_dir:
        out_prot = f"{output_dir}/{output_prefix}.faa"
        out_gff = f"{output_dir}/{output_prefix}.gff"
    else:
        out_prot = f"{output_prefix}.faa"
        out_gff = f"{output_prefix}.gff"
    args = {
        "input": genome,
        "trans_table": str(trans_table),
        "proteins": out_prot,
        "genes": out_gff
    }

    logger.info(" Running prodigal")
    gff = subprocess.run([prodigal, "-i", args["input"], "-g",
                          args["trans_table"], "-m", "-f", "gff",
                          "-a", args["proteins"]],
                         capture_output=True)
    with open(args["genes"], "wb") as f:
        f.write(gff.stdout)
