#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
import argparse
import os


from wiz.qc import qc
# from wiz.misc import util
from wiz.phylogeny import phylogeny

from wiz.version import __version__


def main():
    desc = "take a metagenome-assembled genome and does fancy stuff with it"
    parser = argparse.ArgumentParser(
        prog="wiz",
        description=desc,
        usage="wiz [options] ...",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparsers = parser.add_subparsers(
        title='available subcommands',
        metavar=''
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        default=False,
        help="print software version and exit"
    )

    # phylogeny parser
    parser_phylogeny = subparsers.add_parser(
        'phylogeny',
        prog='wiz phylogeny',
        description='phylogeny module',
        help='wiz: phylogeny module'
    )

    parser_clade = parser_phylogeny.add_mutually_exclusive_group()
    parser_logging_p = parser_phylogeny.add_mutually_exclusive_group()
    parser_logging_p.add_argument(
        '--quiet',
        action='store_true',
        default=False,
        help='Disable info logging'
    )
    parser_logging_p.add_argument(
        '--debug',
        action='store_true',
        default=False,
        help='Enable debug logging'
    )
    parser_phylogeny.add_argument(
        "-g",
        "--genome",
        type=str,
        # metavar="",
        help="input genome bin in fasta format",
        required=True
    )
    parser_clade.add_argument(
        "-c",
        "--clade",
        type=str,
        # metavar="",
        help="clade used for phylogeny and pan-genome analysis"
    )
    parser_clade.add_argument(
        "-d",
        "--clade_dir",
        type=str,
        # metavar="",
        help="do not download a clade and use **.faa file in clade_dir instead"
    )
    parser_phylogeny.add_argument(
        "-o",
        "--output",
        type=str,
        # metavar="",
        default="wiz_output",
        help=f"output directory"
    )
    parser_phylogeny.set_defaults(func=phylogeny.run)

    # qc subparser
    parser_qc = subparsers.add_parser(
        'qc',
        prog='wiz qc',
        description='qc module',
        help='wiz: qc module'
    )
    parser_logging_q = parser_qc.add_mutually_exclusive_group()
    parser_logging_q.add_argument(
        '--quiet',
        action='store_true',
        default=False,
        help='Disable info logging'
    )
    parser_logging_q.add_argument(
        '--debug',
        action='store_true',
        default=False,
        help='Enable debug logging'
    )
    # required arguments
    parser_qc.add_argument(
        "-g",
        "--genomes",
        type=str,
        # metavar="",
        help="input genomes in fasta format",
        required=True,
        nargs='+'
    )
    parser_qc.add_argument(
        "--finchdb",
        type=str,
        required=True,
        help="path to the Finch databank file (.sk)",
        nargs='+'
    )
    parser_qc.add_argument(
        "--taxadb",
        type=str,
        required=True,
        help="Path to the taxadb sqlite database"
    )
    # optional arguments
    parser_qc.add_argument(
        "-o",
        "--output",
        type=str,
        # metavar="",
        default="wiz_output",
        help=f"output directory"
    )
    parser_qc.add_argument(
        "-w",
        "--window",
        type=int,
        default=5000,
        help=f"Size of slicing window for GC statistics"
    )
    parser_qc.add_argument(
        "-c",
        "-cpus",
        type=int,
        default=(os.cpu_count()-1),
        help=f"number of cpus allowed for multi process operations. (default: cpu - 1 = {os.cpu_count()-1})"
    )
    # parser_qc_f = parser_qc.add_mutually_exclusive_group()
    # parser_qc_f.add_argument(
    #     "-a",
    #     "--autofilter",
    #     action='store_true',
    #     default=True,
    #     help=f"Enable auto deletion of contigs with the same ID"
    # )
    # parser_qc_f.add_argument(
    #     "-i",
    #     "--ignorefilter",
    #     action='store_true',
    #     help=f"Disable detection of contigs with the same ID"
    # )
    parser_qc.set_defaults(func=qc.run)

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

        logger = logging.getLogger(__name__)
        logger.info(f"Using wiz version {__version__}")
        logger.debug("Using DEBUG logger")
        args.func(args)
        logging.shutdown()
    except AttributeError as e:
        logger = logging.getLogger(__name__)
        logger.debug(e)
        parser.print_help()
        logging.shutdown()
        #raise  # extra traceback to uncomment for extra debugging powers


if __name__ == "__main__":
    main()
