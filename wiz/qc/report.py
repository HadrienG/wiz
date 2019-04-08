#!/usr/bin/env python
# -*- coding: utf-8

import logging
from math import log
from plotly.offline import plot
from plotly.figure_factory import create_distplot as distplot
from plotly.figure_factory import create_dendrogram as dendrogram
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
    # plot(fig) #* Debug line to plot directly in default internet explorer the graph
    return plot(fig, include_plotlyjs=True, output_type='div')


def scatter_GC_coding_density(bins):
    plotdata = []
    debug = ["name","coding density","GC global"]
    logger.debug(f"{debug[0]:^20}|{debug[1]:^15}|{debug[2]:^15}")
    for item in bins:
        logger.debug(f"{item.name:^20}|{round(item.coding_density,8):>15}|{round(item.gc_global,8):>15}")
        plotdata.append(Scatter(x=[item.coding_density], y=[item.gc_global],name=item.name, mode='markers'))
    layout = Layout(
        title="Coding density on global GC average",
        xaxis=dict(title="coding density (%)", range=[0,100]),
        yaxis=dict(title="Global GC average (%)", range=[0,100])
    )
    fig = Figure(plotdata, layout)
    # plot(fig) #* Debug line to plot directly in default internet explorer the graph
    return plot(fig, include_plotlyjs=True, output_type='div')


def distplot_gc(data):  # waiting a test
    seq_values, seq_names, _, _ = extract_values(data)
    fig = distplot(seq_values, seq_names)
    fig['layout'].update(
        title="Reads ratio per GC average",
        xaxis=dict(title="Average of GC"),
        yaxis=dict(title="Relative amount of reads", range=[0, 1])
    )
    # plot(fig) #* Debug line to plot directly in default internet explorer the graph 
    return plot(fig, include_plotlyjs=True, output_type='div')


def dendrogram_tetra(bins, report):
    id_cells = [Bin.id for Bin in bins]
    cells_value = []
    for id_row in id_cells:
        col_val = []
        for id_col in id_cells:
            if id_col == id_row:
                col_val.append(0)
            else:
                if (id_col, id_row) in report.keys():
                    col_val.append(report[(id_col, id_row)])
                else:
                    col_val.append(report[(id_row, id_col)])
        cells_value.append(col_val)
    cells_value = np.array(cells_value)
    # init fig and create upper dendro
    fig = dendrogram(cells_value, labels=id_cells, orientation='bottom')
    for i in range(len(fig['data'])):    # !
        fig['data'][i]['yaxis'] = 'y2'   # ! found the utility of this sentence
    # create side dendrogram
    fig_side = dendrogram(cells_value, orientation="right")
    for i in range(len(fig_side['data'])):
        fig_side["data"][i]["xaxis"] = "x2"
    # add side dendrogram to fig
    fig.add_traces(fig_side["data"])
    # create heatmap
    dendro_leaves = fig_side['layout']["yaxis"]["ticktext"]
    dendro_leaves = list(map(int, dendro_leaves))
    data_dist = pdist(cells_value)
    heat_data = squareform(data_dist)
    heat_data = heat_data[dendro_leaves, :]
    heat_data = heat_data[:, dendro_leaves]
    heatmap = [
        Heatmap(
            x=dendro_leaves,
            y=dendro_leaves,
            z=heat_data,
            colorscale=[[0, "#080D0F"], [1, "#7BCBF5"]]
            # colorscale=[[0, "#00F532"], [1, "#A80017"]]
        )
    ]
    heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = fig_side['layout']['yaxis']['tickvals']
    # add heatmap data to fig
    fig.add_traces(heatmap)
    # edit layout
    fig['layout'].update({
        #'width': 800,
        #'height': 800,
        'showlegend': False,
        "hovermode": "closest"
        })
    # edit x axis
    fig['layout']['xaxis'].update({
        "domain": [.15, 1],
        "mirror": False,
        "showgrid": False,
        "showline": False,
        "zeroline": False,
        "ticks": ""
    })
    fig['layout'].update({'xaxis2': {
        'domain': [0, .15],
        'mirror': False,
        'showgrid': False,
        'showline': False,
        'zeroline': False,
        'showticklabels': False,
        'ticks': ""
    }})
    # edit yaxis
    fig['layout']["yaxis"].update({
        "domain": [0, .85],
        "mirror": False,
        "showgrid": False,
        "showline": False,
        "zeroline": False,
        "showticklabels": False,
        "ticks": ""
    })
    # edit yaxis2
    fig["layout"].update({'yaxis2': {
        "domain": [.825, .975],
        "mirror": False,
        "showgrid": False,
        "showline": False,
        "zeroline": False,
        "showticklabels": False,
        "ticks": ""
    }})
    fig['layout'].update(
        title="proximity between two sequences based on their composition tetra nucleic"
    )
    # plot(fig)
    return plot(fig, include_plotlyjs=True, output_type='div')
    # ! need to fix xy name of sequence


def contigs_taxonomy(bins):
    figs = []
    dict_link, list_nodes, colors = {}, [], []
    rank_order = [
        "no rank",
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "query"
        ]
    rank_color = {
        "query": "#990000",
        "no rank": "#cc0000",
        "superkingdom": "#ff6600",
        "phylum": "#ff9900",
        "class": "#669900",
        "order": "#3366ff",
        "family": "#330066",
        "genus": "#990099",
        "species": "#cc66cc"
        }
    for genome in bins:
        tax = genome.taxonomy
        for t in tax:
            taxa, jaccard = t
            if taxa["no rank"] != "Not found":
                for rank in rank_order:
                    if rank != "query":
                        if rank in taxa.keys() and taxa[rank] not in list_nodes:
                            list_nodes.append(taxa[rank])
                            colors.append(rank_color[rank])
                    else:
                        if genome.id not in list_nodes:
                            list_nodes.append(genome.id)
                            colors.append(rank_color[rank])
                for rank in rank_order[:-1]:
                    if rank in taxa.keys():
                        actual_node = list_nodes.index(taxa[rank])
                        for r in rank_order[rank_order.index(rank)+1:]:
                            if r in taxa.keys():
                                next_node = list_nodes.index(taxa[r])
                                break
                        if actual_node != next_node:
                            if (actual_node, next_node) in dict_link.keys():
                                dict_link[(actual_node, next_node)] += 1/jaccard
                            else:
                                dict_link[(actual_node, next_node)] = 1/jaccard
                        else:
                            next_node = list_nodes.index(genome.id)
                            if (actual_node, next_node) in dict_link.keys():
                                dict_link[(actual_node, next_node)] += 1/jaccard
                            else:
                                dict_link[(actual_node, next_node)] = 1/jaccard
            # for node in taxa:
            #     if taxa.index(node) == len(taxa)-1:
            #         if (list_nodes.index(node), 0) in dict_link.keys():
            #             dict_link[list_nodes.index(node), 0] += 1/jaccard
            #         else:
            #             dict_link[list_nodes.index(node), 0] = 1/jaccard
            #     else:
            #         next_node = taxa[taxa.index(node)+1]
            #         key = (list_nodes.index(node), list_nodes.index(next_node))
            #         if key in dict_link.keys():
            #             dict_link[key] += 1/jaccard
            #         else:
            #             dict_link[key] = 1/jaccard
    src, target, value = [], [], []
    for k in dict_link:
        s, t = k
        src.append(s)
        target.append(t)
        value.append(dict_link[k])
    data = dict(
        type='sankey',
        arrangement = "freeform",
        node = dict(
            pad = 15,
            #thickness = 15,
            line = dict(
                color = "black",
                width = 0.5
            ),
            label = list_nodes,
            color = colors
        ),
        link = dict(
            source = src,
            target = target,
            value = value
        ))
    layout = dict(
        title = f"Global Taxonomy")
    fig = dict(data=[data], layout=layout)
    plot(fig)
    figs.append(plot(fig, include_plotlyjs=True, output_type='div'))
    return figs


def extract_values(data):
    seq_values, seq_names, bounds, percentil = [], [], [], []
    for dat in data:
        seq_values.append(dat.gc)
        seq_names.append(dat.id)
        bounds.append(dat.gc_bounds)
        percentil.append(dat.gc_percentil)
    return (seq_values, seq_names, bounds, percentil)


def unit(window_size):
    return f"{round_value(window_size)} {factor10_unit(window_size)}"


def factor10_unit(window_size):
    unit = ["b", "Kb", "Mb", "Gb", "Tb"]
    log_value = int(log(window_size, 1000))
    return unit[log_value]


def round_value(window_size):
    log_value = int(log(window_size, 1000))
    numeric_value = round(window_size/1000**log_value, 2)
    return numeric_value


def jinja_report(report_data, args):
    logger.info(" Make a wonderful report for you")
    file_loader = FileSystemLoader("wiz/misc/template", followlinks=True)
    env = Environment(loader=file_loader)
    template = env.get_template("report.html")
    output = template.render(
        param_filein=args.genomes,
        param_fileout=args.output,
        param_cpu = args.c,
        param_window = args.window,
        day_date=time.asctime(),
        average_gc=report_data.gc_scatter_plot,
        gc_dist=report_data.gc_distplot,
        tetra_heatmap=report_data.dendro_tetra)
    return output


def create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        logger.info(" The path folder as been created.")
        logger.info(f" All output files will be stored in {path}")
    else:
        logger.warning(f" The path {path} already exists !")
        logger.warning(" The files in it could have been overwritten!")


def write_QCreport(args, report):
    file_path = os.path.join(args.output, "QC_report.html")
    with open(file_path, "w") as r:
        r.writelines(report)
    logger.info(" The QC report has been successfully written")


class Report:
    def __init__(self, bins, window):
        logger.info(" Develops sophisticated graphs for the report")
        self.gc_scatter_plot = scatter_gc(bins, window)
        self.tetra_distance = distance_calculation(bins)
        self.dendro_tetra = dendrogram_tetra(bins, self.tetra_distance)
        self.gc_distplot = distplot_gc(bins)
        self.gc_coding_density = scatter_GC_coding_density(bins)
        self.contigs_taxonomy = contigs_taxonomy(bins)
