#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
import argparse

from wiz.assets import download


def main():
    desc = "take a metagenome-assembled genome and does fancy stuff with it"
    parser = argparse.ArgumentParser(
        prog="wiz",
        description=desc,
        usage="wiz [options] ...",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
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
    args = parser.parse_args()

    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)

    args.clade = "nostocales"
    a = download.query_assemblies(args.clade, args.output)

    logging.shutdown()


if __name__ == "__main__":
    main()
