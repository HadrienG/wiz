#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import logging

from wiz.misc import util
from wiz.phylogeny import download


def run(args):
    """main function for wiz phylogeny
    """
    logger = logging.getLogger(__name__)

    try:
        # create output dir
        util.create_dir(args.output)

        # download assemblies
        args.clade = "nostocales"
        a = download.query_assemblies(
            args.clade, args.output, quiet=args.quiet, representative=True)

        # create csv report
        output_file = f"{args.output}/assemblies.csv"
        download.create_summary(a, output_file)
    except OSError as e:
        logger.error(f"Directory `{args.output}` already exists. Exiting")
        sys.exit(1)
    except Exception as e:
        raise
    else:
        logger.info("wiz phylogeny finished. Goodbye.")
