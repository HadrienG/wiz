#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
file that gathers the necessary functions for operations
related to drawing graphs.
"""

import logging
from math import log

import numpy
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from plotly.offline import plot
from plotly.figure_factory import create_distplot as distplot
from plotly.figure_factory import create_dendrogram as dendrogram
from plotly.graph_objs import Figure
from plotly.graph_objs import Layout
from plotly.graph_objs import Scatter
from plotly.graph_objs import Heatmap

from wiz.qc.tetra import distance_calculation


logger = logging.getLogger(__name__)


def extract_values(contigs):
    """
    Extracts gc_averages, contig_name, gc_bounds values from contigs
    objects and compiles them in list form to be exploited by chart
    drawing functions
    """
    gc_averages, contig_name, gc_bounds = [], [], []
    for contig in contigs:
        gc_averages.append(contig.gc_cutouts)
        contig_name.append(contig.uid)
        gc_bounds.append(contig.gc_bounds)
    return (gc_averages, contig_name, gc_bounds)


def order_magnitude(number):
    """
    Rounded and assigns the correct order of magnitude to a value.
    """
    return f"{round_value(number)} {magnitude(number)}"


def magnitude(number):
    """
    Set the correct order of magnitude at a given value
    """
    unit = ["b", "Kb", "Mb", "Gb", "Tb", "Pb", "Eb", "Zb", "Yb"]
    log_value = int(log(number, 1000))
    return unit[log_value]


def round_value(number):
    """
    Round a 2-digit numeric value after a comma for an order of magnitude.
    """
    log_value = int(log(number, 1000))
    numeric_value = round(number/1000**log_value, 2)
    return numeric_value


def scatter_gc(contigs, window_size):
    """
    Draw a graph with the proportion of GC along a set of contigs measured
    by x-base segment. The contigs parameter represents a contig object
    list and the window_size parameter represents the size of the clipping
    window that will be displayed at the top of the graph and in the
    abscissa legend. The generated graph is returned in HTML
    """
    logger.debug("drawing scatter_gc graph")

    gc_averages, contigs_name, contigs_bounds = extract_values(contigs)
    plotdata = []

    for averages, name, bounds in zip(
            gc_averages, contigs_name, contigs_bounds):
        # if multiple value in seq
        if isinstance(averages, list):  # type(seq) == list:
            position = [i*window_size for i in range(0, len(averages))]
        # if seq contain a unique value, is not a list but a float
        # we need list for the graph
        else:
            position = [window_size]
            averages = [averages]

        plotdata.append(Scatter(
            x=position,
            y=averages,
            name=name,
            mode='markers',
            marker=dict(size=3)))

        y_down = [bounds[0] for i in range(0, len(averages))]
        y_up = [bounds[1] for i in range(0, len(averages))]

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
        title=f"Average GC per windows of {order_magnitude(window_size)}",
        xaxis=dict(title="Position in the sequence"),
        yaxis=dict(title="Average of GC", range=[0, 100]))

    fig = Figure(plotdata, layout)

    # plot(fig)  #* Debug line to plot directly
    # in default internet explorer the graph

    logger.debug("drawing succesful")

    # return the graph in html format with javascript
    return plot(fig, include_plotlyjs=True, output_type='div')


def scatter_gc_coding_density(contigs):
    """
    Draw a graph representing the average of GC on the proportion of
    coding region of each contig. The generated graph is returned in HTML
    """

    logger.debug("drawing scatter_gc_coding_density graph")

    plotdata = []

    txt_debug = ["name", "coding density", "GC global"]
    logger.debug(f"{txt_debug[0]:^20}|{txt_debug[1]:^15}|{txt_debug[2]:^15}")

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

    # return the graph in html format with javascript
    return plot(fig, include_plotlyjs=True, output_type='div')


def distplot_gc(contigs):
    """
    Draw a graph representing the proportion of contig segments
    with an average of X% GC.
    """

    logger.debug("drawing distplot_gc graph")

    seq_values, seq_names, _ = extract_values(contigs)
    values_temp = []

    for i in seq_values:
        if not isinstance(i, list):
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

    # return the graph in html format with javascript
    return plot(fig, include_plotlyjs=True, output_type='div')


def dendrogram_tetra(contigs, cpu):
    """
    Draws a heatmap representing the proximity of the contigs with respect
    to their tetranucleic composition and sorts them by a hierarchy
    clustering. This graph is generated only if at least 2 contigs are
    submitted. The generated graph is returned in HTML or if it could not
    be generated an information message is returned.
    """

    logger.debug("drawing dendrogram_tetra graph")

    # if more than one contig
    if len(contigs) > 1:
        distances = distance_calculation(contigs, cpu)
        id_cells = [contig.uid for contig in contigs]
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
        cells_value = numpy.array(cells_value)

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

        # edit x axis2
        fig['layout'].update({'xaxis2': {
            'domain': [0, .15],  # side dendrogram proportion in axe x
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': ""
        }})

        # edit y axis
        fig['layout']["yaxis"].update({
            "domain": [0, .85],  # heatmap proportion in axe y
            "mirror": False,
            "showgrid": False,
            "showline": False,
            "zeroline": False,
            "showticklabels": False,
            "ticks": ""
        })

        # edit y axis2
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
        fig['layout'].update(
            title="proximity between two sequences based on their \
                composition tetra nucleic"
        )

        # developement line to draw directly the graph during the run
        # plot(fig)

        logger.debug("drawing succesful")

        # return the graph in html format with javascript
        return plot(fig, include_plotlyjs=True, output_type='div')

    # else :
    logger.debug("drawing unsuccesful, not enough data to draw")
    return "Not enough contigs to draw the heatmap"


def contigs_taxonomy(genome):
    """
    draws a sankey diagram showing the taxonomy of contigs present in
    the bin if these have been identified by taxadb
    """
    logger.debug("drawing contigs_taxonomy graph")

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

    dict_link, nodes, colors = get_nodes(genome, rank_color, rank_order)

    # nodes == false like if len(nodes) != 0
    # if all results are ["no rank"] == "Not found"
    # To fix the problem, check lineage_stepdb or sequence names in sketch file
    if nodes:
        source, target, value = [], [], []

        for keys in dict_link:
            src, targ = keys
            source.append(src)
            target.append(targ)
            value.append(dict_link[keys])

        data = dict(
            type='sankey',
            arrangement="freeform",
            node=dict(
                pad=15,
                line=dict(
                    color="black",
                    width=0.5
                ),
                label=nodes,
                color=colors
            ),
            link=dict(
                source=source,
                target=target,
                value=value
            ))

        layout = dict(
            title=f"Global Taxonomy")

        fig = dict(data=[data], layout=layout)

        # developement line to draw directly the graph during the run
        # plot(fig)

        logger.debug("drawing succesful")

        # return the graph in html format with javascript
        return plot(fig, include_plotlyjs=True, output_type='div')

    return "No taxid found, can be check --lineage_stepdb or --finchdb"


def get_nodes(genome, rank_color, rank_order):
    """
    Function generating the necessary data to draw the taxonomic map.
    """
    dict_link, nodes, colors = {}, [], []
    for contig in genome:
        """
        tax is a tuple list with a dictionary representing the taxonomic
        lineage and jaccard distance of this lineage for the contig
        """
        lineages = contig.taxonomy
        for lineage in lineages:
            lineage_step, jaccard = lineage
            if lineage_step["no rank"] != "Not found":
                for rank in rank_order:
                    if rank != "query":
                        if rank in lineage_step.keys() and \
                                lineage_step[rank] not in nodes:
                            nodes.append(lineage_step[rank])
                            colors.append(rank_color[rank])

                    else:
                        if contig.uid not in nodes:
                            nodes.append(contig.uid)
                            colors.append(rank_color[rank])

                for rank in rank_order[:-1]:
                    if rank in lineage_step.keys():
                        actual_node = nodes.index(lineage_step[rank])

                        for rk in rank_order[rank_order.index(rank)+1:]:
                            if rk in lineage_step.keys():
                                next_node = nodes.index(lineage_step[rk])
                                break

                        if actual_node != next_node:
                            if (actual_node, next_node) in dict_link.keys():
                                dict_link[(actual_node, next_node)] += \
                                    1/jaccard  # More jaccard is small plus
                                # the graph segment must be big, so we
                                # use the inverse of the jaccard value.

                            else:
                                dict_link[(actual_node, next_node)] = 1/jaccard

                        else:
                            next_node = nodes.index(contig.uid)
                            if (actual_node, next_node) in dict_link.keys():
                                dict_link[(actual_node, next_node)] += \
                                    1/jaccard

                            else:
                                dict_link[(actual_node, next_node)] = 1/jaccard
    return dict_link, nodes, colors
