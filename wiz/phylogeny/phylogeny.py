#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import logging

from wiz.misc import path
from wiz.phylogeny import download
from wiz.phylogeny.tools import cd_hit
from wiz.annotate.tools import prodigal


from wiz.misc.path import SoftwareNotFoundError


def run(args):
    """main function for wiz phylogeny
    """
    logger = logging.getLogger(__name__)

    try:
        # create output dir
        path.create_dir(args.output)
        path.create_dir(f"{args.output}/tmp", force=True)

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

        # collect all faa files for cd-hit
        if args.clade:
            protein_files = path.collect(args.output)
        else:
            protein_files = path.collect(
                args.clade_dir) + path.collect(args.output)
            # print(f"faas: {args.clade_dir}")
        # before running cd-hit we need to uncompress input .faa files if they
        # are gzipped
        protein_files = path.uncompress(protein_files)
        print(f"uncompressed:{protein_files}")
        # run cd-hit
        # TODO: cd-hit takes only one file as input.
        # TODO RENAME AND CONCATENATE
        # we'll rename for headers to be of the form >GCFxxxx_PROTID
        cd_hit(protein_files, output_dir=args.output)

    except OSError as e:
        logger.error(f"Directory `{args.output}` already exists. Exiting")
        sys.exit(1)
    except SoftwareNotFoundError as e:
        logger.error(f"`{e}` not in the path. Exiting")
        sys.exit(1)
    else:
        logger.info("wiz phylogeny finished. Goodbye.")
