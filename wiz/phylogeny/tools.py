#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess

from wiz.misc import util


def cd_hit():
    """wrapper function around cd-hit
    """
    # cd-hit -i cyano_proteins.faa -o clustered_proteins_70.faa \
    # -c 0.7 -n 2 -M 100000 -d 0 -T 16
    cd_hit = util.software_exists("cd-hit")
    arguments = {  # placeholders variable
        "input": a,
        "output": b,
        "identity_threshold": c,
        "word_length": d,
        "memory_limit": e,
        "threads": f,
        "-d": 0,
    }
