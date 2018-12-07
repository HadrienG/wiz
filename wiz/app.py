#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
import argparse

from wiz.assets import util
from wiz.assets import download
from wiz.version import __version__


def wiz(args):
    """main function for the wiz software
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Using wiz version {__version__}")
    logger.debug("Using DEBUG logger")

    try:
        # create output dir
        util.create_dir(args.output)

        # download assemblies
        args.clade = "nostocales"
        download_output_dir = f"{args.output}/assemblies"
        util.create_dir(download_output_dir)
        a = download.query_assemblies(
            args.clade, download_output_dir, quiet=args.quiet)

        # create csv report
        output_file = f"{args.output}/assemblies.csv"
        download.create_summary(a, output_file)
    except OSError as e:
        logger.error(f"Directory `{args.output}` exists. Exiting")
        sys.exit(1)
    except Exception as e:
        raise
    else:
        logger.info("wiz finished. Goodbye.")


def main():
    desc = "take a metagenome-assembled genome and does fancy stuff with it"
    parser = argparse.ArgumentParser(
        prog="wiz",
        description=desc,
        usage="wiz [options] ...",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser_logging = parser.add_mutually_exclusive_group()
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        default=False,
        help="print software version and exit"
    )
    parser_logging.add_argument(
        '--quiet',
        action='store_true',
        default=False,
        help='Disable info logging'
    )
    parser_logging.add_argument(
        '--debug',
        action='store_true',
        default=False,
        help='Enable debug logging'
    )
    parser.add_argument(
        "-c",
        "--clade",
        type=str,
        metavar="",
        help="clade used for phylogeny and pan-genome analysis"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        metavar="",
        default="wiz",
        help=f"output directory"
    )
    parser.set_defaults(func=wiz)
    args = parser.parse_args()

    try:
        if args.version:
            print('wiz version %s' % __version__)
            sys.exit(0)
        elif args.quiet:
            logging.basicConfig(level=logging.ERROR)
        elif args.debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        args.func(args)
        logging.shutdown()
    except AttributeError as e:
        logger = logging.getLogger(__name__)
        logger.debug(e)
        parser.print_help()
        logging.shutdown()
        # raise  # extra traceback to uncomment for extra debugging powers


if __name__ == "__main__":
    main()
