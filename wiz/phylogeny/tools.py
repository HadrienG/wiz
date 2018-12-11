#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess

from wiz.misc import path


def cd_hit(proteins, output="clustered_proteins.faa", output_dir=None,
           i=0.40, wl=2):
    """wrapper function around cd-hit
    """
    logger = logging.getLogger(__name__)
    # cd-hit -i cyano_proteins.faa -o clustered_proteins_70.faa \
    # -c 0.7 -n 2 -M 100000 -d 0 -T 16
    cd_hit = path.software_exists("cd-hit")
    args = {  # placeholders variable
        "input": proteins,
        "output": output,
        "identity_threshold": i,
        "word_length": wl,
        "memory_limit": ?,
        "threads": 1,
        "-d": 0,
    }
    logger.info("Running cd-hit")
    subprocess.run([cd_hit, "-i", proteins, "-o", output, "-d", "0", ])
