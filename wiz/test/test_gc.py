#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from wiz.qc import gc


def test_average_gc():  # working in progress
    test_sequences = [
        (Seq("GATCGATGGGCCTATATAGGATCGAAAATC", IUPAC.unambiguous_dna), 46.875),
        (Seq("RATSSATRRSSYTATATARRATYRAAAATY", IUPAC.unambiguous_dna), 9.375)
        ]
    for test in test_sequences:
        test_seq, test_value = test
        assert gc.average_gc(test_seq) == test_value
        assert gc.bio_package_gc(test_seq) == test_value
