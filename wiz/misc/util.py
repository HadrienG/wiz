#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging


def create_dir(dir_name):
    """creates a directory

    Args:
        dir_name (str): name of the directory to create

    Returns:
        str: the absoluta epath of the created directory

    Raises:
        OSError: if directory already exists
    """
    logger = logging.getLogger(__name__)
    try:
        path = os.path.abspath(dir_name)
        os.makedirs(path)
        return path
    except OSError as e:
        logger.debug(f"{dir_name} already exists")
        raise
