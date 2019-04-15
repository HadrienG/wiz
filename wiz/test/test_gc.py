#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from wiz.qc import gc
import pytest
from random import random
from wiz.qc.tools import seq_spliter


def test_average_gc():
    test_sequences = [
        (
            Seq("GATCGATGGGCCTATATAGGATCGAAAATC", IUPAC.unambiguous_dna),
            [43.333333333333336]),
        (
            Seq("RATSSATRRSSYTATATARRATYRAAAATY", IUPAC.unambiguous_dna),
            [13.333333333333334])
        ]
    for test in test_sequences:
        test_seq, test_value = test
        assert gc.average_gc(seq_spliter(
            test_seq, len(test_seq))) == test_value


def test_average_gc_by_frame_with_good_value():
    ch2_Saccharomyces = ""
    ch2_Saccharomyces_result = []
    path_Scerevisae = "wiz/test/saccharomyces_cerevisiae_CH_II.fasta"
    with open(path_Scerevisae, "r") as f:
        _ = f.readline()
        ch2_Saccharomyces = f.read()
        ch2_Saccharomyces.replace("\n", "")
    with open(f"{path_Scerevisae}.result", "r") as result:
        ch2_Saccharomyces_result = [
            entry.replace("\n", "") for entry in result.readlines()]
        ch2_Saccharomyces_result = [
            float(entry) for entry in ch2_Saccharomyces_result]
    test_sequences = [
        (Seq("GATCGATGGGCCTATATAGGATCGAAAATCTAACTTG", IUPAC.unambiguous_dna),
            [60.0, 60.0, 40.0, 40.0, 40.0, 20.0, 20.0], 5),
        (Seq("RATSSATRRSSYTATATARRATYRAAAATY", IUPAC.unambiguous_dna),
            [30.0, 10.0, 0.0], 10)
            ]
    for test in test_sequences:
        if len(test) == 3:
            test_seq, test_value, frame = test
            assert gc.average_gc(seq_spliter(test_seq, frame)) == test_value


def test_average_gc_by_frame_with_bad_value():
    error = [
        "The sequence is void.",
        "The size of the frame is negative or null.",
        "The size of the frame is superior of the sequence length."
    ]

    with pytest.raises(ValueError, message=error[0]):
        gc.average_gc(seq_spliter(Seq("", IUPAC.unambiguous_dna), 10))
    with pytest.raises(ValueError, message=error[1]):
        gc.average_gc(
            seq_spliter(
                Seq(
                    "RATSSATRRSSYTATATARRATYRAAAATY",
                    IUPAC.unambiguous_dna),
                0))
        gc.average_gc(
            seq_spliter(
                Seq(
                    "RATSSATRRSSYTATATARRATYRAAAATY",
                    IUPAC.unambiguous_dna),
                -10))
    with pytest.raises(ValueError, message=error[2]):
        gc.average_gc(
            seq_spliter(
                Seq(
                    "RATSSATRRSSYTATATARRATYRAAAATY",
                    IUPAC.unambiguous_dna),
                5000))


def test_percentil_filter():
    values = [i for i in range(0, 200)]
    short_values = [i for i in range(10, 190)]
    assert gc.percentil_filter(values) == short_values
