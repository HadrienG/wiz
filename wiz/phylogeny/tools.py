#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import subprocess

from wiz.misc import path
from wiz.misc import resources


def cd_hit(proteins, output="clustered_proteins.faa", output_dir=None,
           i=0.40, wl=2, mem=None, cpus=None):
    """wrapper function around cd-hit
    """
    logger = logging.getLogger(__name__)
    # cd-hit -i cyano_proteins.faa -o clustered_proteins_70.faa \
    # -c 0.7 -n 2 -M 100000 -d 0 -T 16

    if mem is None:
        mem = 0.9 * resources.max_mem()
    if cpus is None:
        cpus = resources.max_cpus()

    cd_hit = path.software_exists("cd-hit")
    if output_dir:
        out_prot = f"{output_dir}/{output}"
    else:
        out_prot = output

    input_files = " ".join([str(path) for path in proteins])
    args = {  # placeholders variable
        "input": proteins,
        "output": out_prot,
        "identity": str(i),
        "word_length": str(wl),
        "memory": mem,
        "threads": 1
    }
    logger.info("Running cd-hit")
    logger.warning(f"paths:{proteins}")
    subprocess.run([cd_hit, "-i", input_files, "-o",
                    out_prot, "-d", "0", "-c", args["identity"],
                    "-n", args["word_length"]])
