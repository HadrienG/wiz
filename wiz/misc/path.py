#!/usr/bin/env python
"""contains path utilities for wiz
"""

import os
import gzip
import shutil
import logging

from pathlib import Path


class SoftwareNotFoundError(Exception):
    """Exception to raise when a software is not in the path
    """

    def __init__(self, software):
        super().__init__(f"{software} not found in PATH")


def create_dir(dir_name, force=False):
    """creates a directory

    Args:
        dir_name (str): name of the directory to create
        force (bool): do not raise error if dir exists

    Returns:
        str: the absolute path of the created directory

    Raises:
        OSError: if directory already exists
    """
    logger = logging.getLogger(__name__)
    try:
        path = Path(dir_name)
        path.mkdir(parents=True, exist_ok=force)
        logger.debug(f"created directory: {dir_name}")
        return path.absolute()
    except OSError as e:
        logger.debug(f"{dir_name} already exists")
        raise


def software_exists(software_name):
    """check if a command-line utility exists and is in the path
    """
    s = shutil.which(software_name)
    if s is not None:
        return s
    else:
        raise SoftwareNotFoundError(software_name)


def collect(dir_name, extension=".faa.gz", recursive=False):
    """collect all files of the same extension in a directory
    """
    extension_uncompressed = extension.split(".")[:-1]
    path = Path(dir_name)
    result_compressed = path.glob(f"**/*{extension}")
    result_uncompressed = path.glob(f"**/*{extension_uncompressed}")
    return list(result_compressed) + list(result_uncompressed)


def uncompress(paths, remove=True):
    """uncompress all .gz files from a list of paths
    """
    uncompressed_paths = []
    for path in paths:
        out_path = path.with_suffix("")
        if is_gzipped(path):
            with gzip.open(path, "rb") as i, open(out_path, "wb") as o:
                shutil.copyfileobj(i, o)
            if remove:
                Path.unlink(path)
            uncompressed_paths.append(out_path)
        else:
            uncompressed_paths.append(path)
    return(uncompressed_paths)


def is_gzipped(infile):
    """Check in a file is in gzip format or not
    """
    logger = logging.getLogger(__name__)

    magic_number = b'\x1f\x8b'
    f = open(infile, 'rb')
    with f:
        try:
            assert f.read(2) == magic_number
        except AssertionError as e:
            logger.debug(f'{infile} is not gzipped')
            return False
        else:
            logger.debug(f'{infile} is gzipped')
            return True
