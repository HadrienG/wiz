#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil

from wiz.misc import util
from wiz.phylogeny import download


def setup_module():
    util.create_dir("test_output_dir")


def teardown_module():
    shutil.rmtree("test_output_dir")


def test_query_assemblies():
    assemblies = download.query_assemblies(
        "Nostoc", "test_output_dir", representative=True)
    # bad test, will break if a new genome is added
    assert len(assemblies) == 3
