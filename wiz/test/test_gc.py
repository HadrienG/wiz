#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from wiz.qc import gc
import pytest
from random import random


def test_average_gc():
    test_sequences = [
        (Seq("GATCGATGGGCCTATATAGGATCGAAAATC", IUPAC.unambiguous_dna), 46.875),
        (Seq("RATSSATRRSSYTATATARRATYRAAAATY", IUPAC.unambiguous_dna), 9.375)
        ]
    for test in test_sequences:
        test_seq, test_value = test
        assert gc.average_gc(test_seq) == test_value
        assert gc.bio_package_gc(test_seq) == test_value


def test_average_gc_by_frame_with_good_value():
    ch2_Saccharomyces = ""
    ch2_Saccharomyces_result = []
    with open("saccharomyces cerevisiae CH II.fasta", "r") as f:
        _ = f.readline()
        ch2_Saccharomyces = f.read()
        ch2_Saccharomyces.replace("\n", "")
    with open("saccharomyces cerevisiae CH II.fasta.result", "r") as result:
        ch2_Saccharomyces_result = [entry.replace("\n", "") for entry in result.readlines()]
        ch2_Saccharomyces_result = [float(entry) for entry in ch2_Saccharomyces_result]
    test_sequences = [
        (Seq("GATCGATGGGCCTATATAGGATCGAAAATCTAACTTG", IUPAC.unambiguous_dna),
            [(60, "GATCG"), (60, "ATGGG"), (40, "CCTAT"), (40, "ATAGG"), (40, "ATCGA"), (20, "AAATC"), (20, "TAACT"),
                (50, "TG")], 5),
        (Seq("RATSSATRRSSYTATATARRATYRAAAATY", IUPAC.unambiguous_dna),
            [(30, "RATSSATRRS"), (10, "SYTATATARR"), (0, "ATYRAAAATY")], 10)
            ]
    for test in test_sequences:
        if len(test) == 3:
            test_seq, test_value, frame = test
            assert gc.average_gc_by_frame(test_seq, frame) == test_value


def test_average_gc_by_frame_with_bad_value():
    with pytest.raises(ValueError, message="The sequence is void"):
        gc.average_gc_by_frame(Seq("", IUPAC.unambiguous_dna))
    with pytest.raises(ValueError, message="The size of the frame is negative or null."):
        gc.average_gc_by_frame(Seq("RATSSATRRSSYTATATARRATYRAAAATY", IUPAC.unambiguous_dna), 0)
        gc.average_gc_by_frame(Seq("RATSSATRRSSYTATATARRATYRAAAATY", IUPAC.unambiguous_dna), -10)
    with pytest.raises(ValueError, message="The size of the frame is superior of the sequence length."):
        gc.average_gc_by_frame(Seq("RATSSATRRSSYTATATARRATYRAAAATY", IUPAC.unambiguous_dna))


def test_deviating_value():
    test_list = [random()*100 for i in range(0, S250)]
    test1 = gc.deviating_contigs(test_list)
    assert len(test1) == 13
    test2 = gc.deviating_contigs(test_list, 80)
    assert len(test2) == 50
    for value in test1:
        assert value in test2
