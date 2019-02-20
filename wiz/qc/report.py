#!/usr/bin/env python
# -*- coding: utf-8

import logging
from math import log
from plotly.offline import plot
from plotly.figure_factory import create_distplot as distplot
from plotly.graph_objs import Histogram, Figure, Layout, Scatter

logger = logging.getLogger(__name__)


class Report:
    def __init__(self, bins, window):
        self.gc_scatter_plot = scatter_gc(bins, window)
        self.gc_distplot = ""  # distplot_gc(bins)
        self.tetra_distribution = ""

    # def __repr__(self):
    #     return self.gc_scatter_plot + self.gc_distplot


def scatter_gc(data, window_size):  # waiting a test
    seq_values, seq_names, seq_bounds, seq_percentil = extract_values(data)
    plotdata = []
    for seq, name, bound, percentil in zip(seq_values, seq_names, seq_bounds, seq_percentil):
        position = [i*window_size for i in range(0, len(seq))]
        plotdata.append(Scatter(x=position, y=seq, name=name, mode='markers', marker=dict(size=3)))
        y_down = [bound[0] for i in range(0, len(seq))]
        y_up = [bound[1] for i in range(0, len(seq))]
        plotdata.append(Scatter(x=position, y=y_down, name=name+" "+str(percentil[0])+"e percentil", mode='lines', line=dict(width=1)))
        plotdata.append(Scatter(x=position, y=y_up, name=name+" "+str(percentil[1])+"e percentil", mode='lines', line=dict(width=1)))
    layout = Layout(  # * Try to change scatter in plot or bar
        title=f"Average GC per windows of {unit(window_size)}",
        xaxis=dict(title="Position in the sequence"),
        yaxis=dict(title="Average of GC", range=[0, 100]))
    fig = Figure(plotdata, layout)
    plot(fig)
    # return plot(fig, include_plotlyjs=True, output_type='div')


def distplot_gc(data):  # waiting a test
    seq_values, seq_names, _,_ = extract_values(data)
    fig = distplot(seq_values, seq_names)
    fig['layout'].update(
        title="Reads ratio per GC average",
        xaxis=dict(title="Average of GC"),
        yaxis=dict(title="Relative amount of reads", range=[0, 1])
    )
    # return plot(fig, include_plotlyjs=True, output_type='div')
    plot(fig)  # just here to help in the dev of this function
# TODO comment the displot graph


def extract_values(data):  # I think it's OK
    seq_values, seq_names, bounds, percentil = [], [], [], []
    for dat in data:
        seq_values.append(dat.gc)
        seq_names.append(dat.id)
        bounds.append(dat.gc_bounds)
        percentil.append(dat.gc_percentil)
    return (seq_values, seq_names, bounds, percentil)


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
