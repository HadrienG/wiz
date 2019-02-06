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


def test_average_gc_by_frame():
    ch2_Saccharomyces = ""
    ch2_Saccharomyces_result = []
    with open("test_gc.fasta", "r") as f:
        _ = f.readline()
        ch2_Saccharomyces = f.read()
        ch2_Saccharomyces.replace("\n", "")
    with open("test_gc.result", "r") as result:
        ch2_Saccharomyces_result = [entry.replace("\n", "") for entry in result.readlines()]
        ch2_Saccharomyces_result = [float(entry) for entry in ch2_Saccharomyces_result]
    #   ch2_Saccharomyces_result = [float(entry.replace("\n","")) for entry in result.readlines()]   mistakes here
    test_sequences = [
        (Seq("GATCGATGGGCCTATATAGGATCGAAAATCTAACTTG", IUPAC.unambiguous_dna),
            [60, 60, 40, 40, 40, 20, 20,
                50], 5),
        (Seq("RATSSATRRSSYTATATARRATYRAAAATY", IUPAC.unambiguous_dna),
            [30, 10, 0], 10),
        (Seq(ch2_Saccharomyces, IUPAC.unambiguous_dna), ch2_Saccharomyces_result)  # ,
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
