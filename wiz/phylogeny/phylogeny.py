#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import logging

from wiz.misc import path
from wiz.phylogeny import download
from wiz.annotate.tools import prodigal

from wiz.misc.path import SoftwareNotFoundError


def run(args):
    """main function for wiz phylogeny
    """
    logger = logging.getLogger(__name__)

    try:
        # create output dir
        path.create_dir(args.output)

        # download assemblies
        if args.clade:
            args.clade = "nostocales"
            a = download.query_assemblies(
                args.clade, args.output, quiet=args.quiet, representative=True)

            # create csv report
            output_file = f"{args.output}/assemblies.csv"
            download.create_summary(a, output_file)

        elif not args.clade_dir:
            logger.error("One of --clade or --clade_dir is required. Exiting")
            sys.exit(1)

        # run prodigal
        prodigal(args.genome, output_dir=args.output)
        if args.clade:
            protein_files = path.collect(args.output)
        else:
            protein_files = path.collect(
                args.clade_dir) + path.collect(args.output)

        # run cd-hit
        # ...
    except OSError as e:
        logger.error(f"Directory `{args.output}` already exists. Exiting")
        sys.exit(1)
    except SoftwareNotFoundError as e:
        logger.error(f"`{e}` not in the path. Exiting")
        sys.exit(1)
    else:
        logger.info("wiz phylogeny finished. Goodbye.")
