#!/usr/bin/env python
# -*- coding: utf-8

import logging
from math import log
from plotly.offline import plot
from plotly.figure_factory import create_distplot as distplot
from plotly import figure_factory as ff
from plotly.graph_objs import Histogram, Figure, Layout, Scatter, Heatmap
from wiz.qc.tetra import distance_calculation
from jinja2 import Environment, FileSystemLoader
import os
import time
import numpy as np
from scipy.spatial.distance import pdist, squareform

logger = logging.getLogger(__name__)


def scatter_gc(data, window_size):  # waiting a test
    seq_values, seq_names, seq_bounds, seq_percentil = extract_values(data)
    plotdata = []
    for seq, name, bound, percentil in zip(
            seq_values, seq_names, seq_bounds, seq_percentil):
        position = [i*window_size for i in range(0, len(seq))]
        plotdata.append(Scatter(x=position, y=seq, name=name,
            mode='markers', marker=dict(size=3)))
        y_down = [bound[0] for i in range(0, len(seq))]
        y_up = [bound[1] for i in range(0, len(seq))]
        plotdata.append(Scatter(x=position, y=y_down,
            name=name+" "+str(percentil[0])+"e percentil",
            mode='lines', line=dict(width=1), opacity=0.25, showlegend=False))
        plotdata.append(Scatter(x=position, y=y_up,
            name=name+" "+str(percentil[1])+"e percentil", mode='lines',
            line=dict(width=1), opacity=0.25, showlegend=False))
    layout = Layout(  # * Try to change scatter in plot or bar
        title=f"Average GC per windows of {unit(window_size)}",
        xaxis=dict(title="Position in the sequence"),
        yaxis=dict(title="Average of GC", range=[0, 100]))
    fig = Figure(plotdata, layout)
    #plot(fig)
    return plot(fig, include_plotlyjs=True, output_type='div')


def distplot_gc(data):  # waiting a test
    seq_values, seq_names, _, _ = extract_values(data)
    fig = distplot(seq_values, seq_names)
    fig['layout'].update(
        title="Reads ratio per GC average",
        xaxis=dict(title="Average of GC"),
        yaxis=dict(title="Relative amount of reads", range=[0, 1])
    )
    # plot(fig)  # just here to help in the dev of this function
    return plot(fig, include_plotlyjs=True, output_type='div')
# TODO comment the displot graph


def dendro_tetra(bins,report):
    id_cell = [bin.id for bin in bins]
    cell_values = []
    for id_row in id_cell:
        col_value = []
        for id_col in id_cell:
            if id_col == id_row:
                col_value.append(0)
            elif (id_col, id_row) in report.keys():
                col_value.append(report[(id_col, id_row)])
            else:
                col_value.append(report[(id_row, id_col)])
        cell_values.append(col_value)
    data_array = np.array(cell_values)
    labels = id_cell
    #src : https://plot.ly/python/dendrogram/
    # Initialize figure by creating upper dendrogram
    figure = ff.create_dendrogram(data_array, orientation='bottom', labels=labels)
    for i in range(len(figure['data'])):
        figure['data'][i]['yaxis'] = 'y2'

    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(data_array, orientation='right')
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'
    figure['data']=list(figure['data'])
    # Add Side Dendrogram Data to Figure
    figure['data'].extend(dendro_side['data'])

    # Create Heatmap
    dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))
    data_dist = pdist(data_array)
    heat_data = squareform(data_dist)
    heat_data = heat_data[dendro_leaves,:]
    heat_data = heat_data[:,dendro_leaves]

    heatmap = [
        go.Heatmap(
            x = dendro_leaves,
            y = dendro_leaves,
            z = heat_data,
            colorscale = 'YIGnBu'
        )
    ]

    heatmap[0]['x'] = figure['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    figure['data'].extend(heatmap)

    # Edit Layout
    figure['layout'].update({'width':800, 'height':800,
                            'showlegend':False, 'hovermode': 'closest',
                            })
    plot(fig)


def table_tetra(bins, report):
    id_cell = [bin.id for bin in bins]
    cell_values = []
    for id_row in id_cell:
        col_value = []
        for id_col in id_cell:
            if id_col == id_row:
                col_value.append(0)
            elif (id_col, id_row) in report.keys():
                col_value.append(report[(id_col, id_row)])
            else:
                col_value.append(report[(id_row, id_col)])
        cell_values.append(col_value)
    trace = Heatmap(
        z=cell_values,
        x=id_cell,
        y=id_cell,
        # colorscale=[[0,"#0CFF15"],[0.5,"#FFC900"],[1,"#700F00"]]
        colorscale=[[0, "#022A33"], [1, "#7BCBF5"]]
    )
    data = [trace]
    # plot(data)  # dev line
    return plot(data, include_plotlyjs=True, output_type='div')


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


class Report:
    def __init__(self, bins, window):
        logger.info("Make a wonderful report for you")
        self.gc_scatter_plot = scatter_gc(bins, window)
        self.gc_distplot = distplot_gc(bins)
        self.tetra_distance = distance_calculation(bins)
        self.tetra_heatmap = table_tetra(bins, self.tetra_distance)
        dendro_tetra(bins, self.tetra_distance)

    # def __repr__(self):
    #     return self.gc_scatter_plot + self.gc_distplot


def jinja_report(report_data):
    file_loader = FileSystemLoader("wiz/misc/template", followlinks=True)
    env = Environment(loader=file_loader)
    template = env.get_template("report.html")
    output = template.render(
        day_date=time.asctime(),
        average_gc=report_data.gc_scatter_plot,
        gc_dist=report_data.gc_distplot,
        tetra_heatmap=report_data.tetra_heatmap)
    return output


def create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        logger.info(" The path folder as been created.")
        logger.info(f" All output files will be stored in {path}")
    else:
        logger.warning(f" The path {path} already exists !")
        logger.warning(" The files in it could have been overwritten!")


def write_QCreport(path, report):
    file_path = os.path.join(path, "QC_report.html")
    with open(file_path, "w") as r:
        r.writelines(report)
    logger.info(" The QC report has been successfully written")
