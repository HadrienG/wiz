#!/usr/bin/env python
# -*- coding: utf-8

import logging
from math import log
from plotly.offline import plot
from plotly.figure_factory import create_distplot as distplot
from plotly.graph_objs import Histogram, Figure, Layout, Scatter

logger = logging.getLogger(__name__)


def scatter_gc(data, window_size=5000):  # waiting a test
    seq_values, seq_names, bounds = extract_values(data)
    plotdata = []
    for seq, name,bound in zip(seq_values, seq_names, bounds):
        position = [i*window_size for i in range(0, len(seq))]
        plotdata.append(Scatter(x=position, y=seq, name=name, mode='markers'))
        y_down = [bound[0] for i in range(0, len(seq))]
        y_up = [bound[1] for i in range(0, len(seq))]
        plotdata.append(Scatter(x=position, y=y_down, name=name+" bounds_DOWN", mode='markers'))
        plotdata.append(Scatter(x=position, y=y_up, name=name+" bounds_UP", mode='markers'))
    for bound in bounds:
        position = [i*window_size for i in range(0, len(seq))]
        
    layout = Layout(  # * Try to change scatter in plot or bar
        title=f"Average GC per windows of {unit(window_size)}",
        xaxis=dict(title=f"Pourcent of GC"),
        yaxis=dict(title="Sequence number", range=[0, 100]))
    fig = Figure(plotdata, layout)
    plot(fig)


def distplot_gc(data):  # waiting a test
    seq_values, seq_names = extract_values(data)
    fig = distplot(seq_values, seq_names)
    plot(fig)
# TODO comment the displot graph


def extract_values(data):  # I think it's OK
    seq_values, seq_names, bounds = [], [], []
    for dat in data:
        seq_values.append(dat.gc)
        seq_names.append(dat.id)
        bounds.append(dat.gc_bounds)
    return (seq_values, seq_names, bounds)


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
