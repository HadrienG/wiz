#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from wiz.qc import gc


def test_average_gc():
    test_sequences = [
        (Seq("GATCGATGGGCCTATATAGGATCGAAAATC", IUPAC.unambiguous_dna), 46.875),
        (Seq("RATSSATRRSSYTATATARRATYRAAAATY", IUPAC.unambiguous_dna), 9.375)
        ]
    for test in test_sequences:
        test_seq, test_value = test
        assert gc.average_gc(test_seq) == test_value
        assert gc.bio_package_gc(test_seq) == test_value


def test_average_gc_by_frame:
    test_sequences = [
        (Seq("GATCGATGGGCCTATATAGGATCGAAAATCTAACTTG", IUPAC.unambiguous_dna),
            ["GATCG", "ATGGG", "CCTAT", "ATAGG", "ATCGA", "AAATC", "TAACT",
                "TG"]),
        (Seq("RATSSATRRSSYTATATARRATYRAAAATY", IUPAC.unambiguous_dna),
            ["RATSSATRRS", "SYTATATARR", "ATYRAAAATY"], 10)  # ,
            #  (Seq("RATSSATRRSSYTATATARRATYRAAAATY", IUPAC.unambiguous_dna), ValueError, 0),
            #  (Seq("RATSSATRRSSYTATATARRATYRAAAATY", IUPAC.unambiguous_dna), ValueError, -10),
            #  (Seq("", IUPAC.unambiguous_dna), ValueError, 10)
            ]
    for test in test_sequences:
        if len(test) == 3:
            test_seq, test_value, frame = test
            assert gc.average_gc_by_frame(test_seq, frame) == test_value
        else:
            test_seq, test_value = test
            assert gc.average_gc_by_frame(test_seq) == test_value
        