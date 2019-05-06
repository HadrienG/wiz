#!/usr/bin/env python3
# -*- coding: utf-8

"""
Apply hmmpress to all files passed in parameter with -f .
"""

import shutil
import subprocess
import argparse


def hmmpress(files):
    hmmpress = software_exists("hmmpress")
    for file in files:
        s_args = [hmmpress, "-f", file]
        subprocess.run(s_args)


def software_exists(software_name):
    """check if a command-line utility exists and is in the path
    """ # this portion of code was wrote by Hadrien Gourle
    s = shutil.which(software_name)
    if s is not None:
        return s
    else:
        raise SoftwareNotFoundError(software_name)


def main(args):
    hmmpress(args.files)
    print("done")

parser = argparse.ArgumentParser(
    prog="HMMPressMulti",
    usage="HMMPM -f file1.hmm, file2.hmm, ...",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
    '-f',
    "--files",
    type=str,
    metavar=" [str]",
    help="input files to hmmpress",
    required=True,
    nargs="+"
)

args = parser.parse_args()

if __name__ == "__main__":
    main(args)
