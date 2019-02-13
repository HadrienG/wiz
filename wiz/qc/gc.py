#!/usr/bin/env python
# -*- coding:utf-8 -*-

from math import log
from Bio.SeqUtils import GC
from plotly.offline import plot
from plotly.figure_factory import create_distplot as distplot
from plotly.graph_objs import Histogram, Figure, Layout, Scatter
from numpy import percentile


def check_window_size(sequence, window_size):       # I think it's OK
    if len(sequence) == 0:
        error = "The sequence is void."
        raise ValueError(error)
    if window_size < 0:
        error = "The size of the window is negative."
        raise ValueError(error)
    if window_size == 0:
        error = "The size of the window is null."
        raise ValueError(error)
    if window_size > len(sequence):
        error = "The size of the window is superior of the sequence length."
        raise ValueError(error)
    return True


def average_gc(sequence, window_size=5000, truncate=False):  # I think it's OK
    """
    Return a list with for each window of size window_size the average of GC.
    The last one windows can be more small that the others windows due of the
    sequence length. The result of the window of length inferior at window_size
    can be removed with the parameter truncate=True
    """
    if check_window_size(sequence, window_size):
        average = []
        seq_size = len(sequence)
        for pos in range(0, seq_size, window_size):
            if pos+window_size < seq_size:
                subseq = sequence[pos:pos+window_size]
            else:
                if not truncate:
                    subseq = sequence[pos:]
                else:
                    break
            average.append(GC(subseq))
        return average


def cleaning(average):  # WIP need testing
    clean_seq = percentile(average, 0.95)
    return clean_seq


# Graph functions

def distplot_gc(data):  # waiting a test
    seq_values, seq_names = extract_values(data)
    fig = distplot(seq_values, seq_names)
    plot(fig)


def histogram_gc(data, window_size=5000):  # waiting a test
    seq_values, seq_names = extract_values(data)
    hist_data = []
    for values, name in zip(seq_values, seq_names):
        hist_data.append(Histogram(
            x=values,
            name=name,
            opacity=0.5
        ))
    layout = Layout(
        title=f"Average GC per windows of {unit(window_size)}",
        xaxis=dict(title=f"Pourcent of GC"),
        yaxis=dict(title="Sequence number"),
        barmode='overlay'
    )
    fig = Figure(data=hist_data, layout=layout)
    plot(fig)


def scatter_gc(data, window_size=5000):  # waiting a test
    seq_values, seq_names = extract_values(data)
    plotdata = [Scatter(x=[0], y=[0], name="bound down"), Scatter(x=[0], y=[100], name="bound up")]
    for seq, name in zip(seq_values, seq_names):
        position = [i*window_size for i in range(0, len(seq))]
        plotdata.append(Scatter(x=position, y=seq, name=name))
    plot(plotdata)

# facultative functions =============================================


def unit(window_size):  # I think it's OK
    return f"{round_value(window_size)} {factor10_unit(window_size)}"


def factor10_unit(window_size):  # I think it's OK
    unit = ["b", "Kb", "Mb", "Gb", "Tb"]
    log_value = int(log(window_size, 1000))
    return unit[log_value]


def round_value(window_size):  # I think it's OK
    log_value = int(log(window_size, 1000))
    numeric_value = round(window_size/1000**log_value, 2)
    return numeric_value


def extract_values(data):  # I think it's OK
    seq_values, seq_names = [], []
    for dat in data:
        seq, name = dat
        seq_values.append(seq)
        seq_names.append(name)
    return (seq_values, seq_names)
