#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from plotly.offline import plot
from plotly.figure_factory import create_distplot as distplot
from plotly.figure_factory import create_dendrogram as dendrogram
from plotly.graph_objs import Histogram, Figure, Layout, Scatter, Heatmap
from wiz.qc.tetra import distance_calculation
import numpy as np
from scipy.spatial.distance import pdist, squareform
from math import log

logger = logging.getLogger(__name__)


def extract_values(contigs):
    seq_values, seq_names, bounds = [], [], []
    for contig in contigs:
        seq_values.append(contig.gc_cutouts)
        seq_names.append(contig.id)
        bounds.append(contig.gc_bounds)
    return (seq_values, seq_names, bounds)


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


def scatter_gc(contigs, window_size):
    logger.debug("drawing scatter_gc graph")
    seq_values, seq_names, seq_bounds = extract_values(contigs)
    plotdata = []
    for seq, name, bound in zip(
            seq_values, seq_names, seq_bounds):
        # if multiple value in seq
        if type(seq) == list:
            position = [i*window_size for i in range(0, len(seq))]
        # if seq contain a unique value, is not a list bat a float
        # we need list for the graph
        else:
            position = [window_size]
            seq = [seq]
        plotdata.append(Scatter(
            x=position,
            y=seq,
            name=name,
            mode='markers',
            marker=dict(size=3)))
        y_down = [bound[0] for i in range(0, len(seq))]
        y_up = [bound[1] for i in range(0, len(seq))]
        plotdata.append(Scatter(
            x=position,
            y=y_down,
            name=name+" bounds",
            mode='lines',
            line=dict(width=1),
            opacity=0.25,
            showlegend=False))
        plotdata.append(Scatter(
            x=position,
            y=y_up,
            name=name+" bounds",
            mode='lines',
            line=dict(width=1),
            opacity=0.25,
            showlegend=False))
    layout = Layout(  # * Try to change scatter in plot or bar
        title=f"Average GC per windows of {unit(window_size)}",
        xaxis=dict(title="Position in the sequence"),
        yaxis=dict(title="Average of GC", range=[0, 100]))
    fig = Figure(plotdata, layout)
    # plot(fig)  #* Debug line to plot directly
    # in default internet explorer the graph
    logger.debug("drawing succesful")
    return plot(fig, include_plotlyjs=True, output_type='div')


def scatter_GC_coding_density(contigs):
    logger.debug("drawing scatter_gc_coding_density graph")
    plotdata = []
    debug = ["name", "coding density", "GC global"]
    logger.debug(f"{debug[0]:^20}|{debug[1]:^15}|{debug[2]:^15}")
    for contig in contigs:
        txt = (contig.name, contig.coding_density, contig.gc_contig)
        logger.debug(
            f"{txt[0]:^20}|{round(txt[1], 8):>15}|{round(txt[2], 8):>15}")
        plotdata.append(Scatter(
            x=[contig.coding_density],
            y=[contig.gc_contig],
            name=contig.name,
            mode='markers'))
    layout = Layout(
        title="Coding density on global GC average",
        xaxis=dict(title="coding density (%)", range=[0, 100]),
        yaxis=dict(title="Global GC average (%)", range=[0, 100])
    )
    fig = Figure(plotdata, layout)
    # plot(fig) #* Debug line to plot directly
    # in default internet explorer the graph
    logger.debug("drawing succesful")
    return plot(fig, include_plotlyjs=True, output_type='div')


def distplot_gc(contigs):
    logger.debug("drawing distplot_gc graph")
    seq_values, seq_names, _ = extract_values(contigs)
    values_temp = []
    for i in seq_values:
        if type(i) != list:
            i = [i, 1]
        values_temp.append(i)
    fig = distplot(values_temp, seq_names, show_hist=False)
    fig['layout'].update(
        title="Reads ratio per GC average",
        xaxis=dict(title="Average of GC", range=[20, 80]),
        yaxis=dict(title="Relative amount of reads", range=[0, 1])
    )
    # plot(fig) #* Debug line to plot directly
    # in default internet explorer the graph
    logger.debug("drawing succesful")
    return plot(fig, include_plotlyjs=True, output_type='div')


def dendrogram_tetra(contigs, cpu):
    logger.debug("drawing dendrogram_tetra graph")
    # if more than one contig
    if len(contigs) > 1:
        distances = distance_calculation(contigs, cpu)
        id_cells = [contig.id for contig in contigs]
        cells_value = []
        for id_row in id_cells:
            col_val = []
            for id_col in id_cells:
                if id_col == id_row:
                    col_val.append(0)
                else:
                    if (id_col, id_row) in distances.keys():
                        col_val.append(distances[(id_col, id_row)])
                    else:
                        col_val.append(distances[(id_row, id_col)])
            cells_value.append(col_val)
        cells_value = np.array(cells_value)
        # init fig and create upper dendro
        fig = dendrogram(cells_value, labels=id_cells, orientation='bottom')
        # create the colorscale legend
        for i in range(len(fig['data'])):
            fig['data'][i]['yaxis'] = 'y2'
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
            # 'width': 800,
            # 'height': 800,
            'showlegend': False,
            "hovermode": "closest"
            })
        # edit x axis
        fig['layout']['xaxis'].update({
            "domain": [.15, 1],  # heatmap proportion in axe x
            "mirror": False,
            "showgrid": False,
            "showline": False,
            "zeroline": False,
            "ticks": ""
        })
        fig['layout'].update({'xaxis2': {
            'domain': [0, .15],  # side dendrogram proportion in axe x
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        }})
        # edit yaxis
        fig['layout']["yaxis"].update({
            "domain": [0, .85],  # heatmap proportion in axe y
            "mirror": False,
            "showgrid": False,
            "showline": False,
            "zeroline": False,
            "showticklabels": False,
            "ticks": ""
        })
        # edit yaxis2
        fig["layout"].update({'yaxis2': {
            "domain": [.825, .975],  # upper dendrogram proportion in axe x
            "mirror": False,
            "showgrid": False,
            "showline": False,
            "zeroline": False,
            "showticklabels": False,
            "ticks": ""
        }})
        # edit title of the graph
        title = "proximity between two sequences based "
        title += "on their composition tetra nucleic"
        fig['layout'].update(
            title=title
        )
        # developement line to draw directly the graph during the run
        # plot(fig)

        # draw the graph in html format with javascript
        return plot(fig, include_plotlyjs=True, output_type='div')
        logger.debug("drawing succesful")
    else:
        logger.debug("drawing unsuccesful, not enough data to draw")
        return("Not enough data to draw the map")


def contigs_taxonomy(bins):
    logger.debug("drawing contigs_taxonomy graph")
    dict_link, nodes, colors = {}, [], []
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
                        if rank in taxa.keys() and taxa[rank] not in nodes:
                            nodes.append(taxa[rank])
                            colors.append(rank_color[rank])
                    else:
                        if genome.id not in nodes:
                            nodes.append(genome.id)
                            colors.append(rank_color[rank])
                for rank in rank_order[:-1]:
                    if rank in taxa.keys():
                        actual_node = nodes.index(taxa[rank])
                        for r in rank_order[rank_order.index(rank)+1:]:
                            if r in taxa.keys():
                                next_node = nodes.index(taxa[r])
                                break
                        if actual_node != next_node:
                            if (actual_node, next_node) in dict_link.keys():
                                dict_link[(
                                    actual_node,
                                    next_node)] += 1/jaccard
                            else:
                                dict_link[(actual_node, next_node)] = 1/jaccard
                        else:
                            next_node = nodes.index(genome.id)
                            if (actual_node, next_node) in dict_link.keys():
                                dict_link[(
                                    actual_node,
                                    next_node)] += 1/jaccard
                            else:
                                dict_link[(actual_node, next_node)] = 1/jaccard
    if len(nodes) != 0:
        src, target, value = [], [], []
        for k in dict_link:
            s, t = k
            src.append(s)
            target.append(t)
            value.append(dict_link[k])
        data = dict(
            type='sankey',
            arrangement="freeform",
            node=dict(
                pad=15,
                # thickness = 15,
                line=dict(
                    color="black",
                    width=0.5
                ),
                label=nodes,
                color=colors
            ),
            link=dict(
                source=src,
                target=target,
                value=value
            ))
        layout = dict(
            title=f"Global Taxonomy")
        fig = dict(data=[data], layout=layout)
        # plot(fig)
        logger.debug("drawing succesful")
        return plot(fig, include_plotlyjs=True, output_type='div')
    else:
        return "No taxid found, can be check --taxadb or --finchdb"
