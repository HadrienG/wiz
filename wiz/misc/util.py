#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
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

    Returns:
        str: the absolute path of the created directory

    Raises:
        OSError: if directory already exists
    """
    logger = logging.getLogger(__name__)
    try:
        path = Path(dir_name)
        path.mkdir(exist_ok=force)

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
