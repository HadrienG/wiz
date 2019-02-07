#!/usr/bin/env python
# -*- coding:utf-8 -*-

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from plotly.offline import plot
import plotly.figure_factory as ff
import plotly.graph_objs as go
from numpy import histogram


def check_frame_size(sequence, len_frame):
    if len(sequence) == 0:
        error = "The sequence is void."
        raise ValueError(error)
    if len_frame <= 0:
        error = "The size of the frame is negative or null."
        raise ValueError(error)
    if len_frame > len(sequence):
        error = "The size of the frame is superior of the sequence length."
        raise ValueError(error)
    return True


def average_gc_by_frame(sequence, len_frame=5000):
    if check_frame_size(sequence, len_frame):
        average_gc_by_windows = []
        for pos in range(0, len(sequence), len_frame):
            if pos+len_frame < len(sequence):
                subsequence = sequence[pos:pos+len_frame]
            else:
                subsequence = sequence[pos:]
            average_gc = GC(subsequence)
            average_gc_by_windows.append(average_gc, subsequence)
        return average_gc_by_windows


def gc_hist(sequence_list, len_frame=5000):
    """
    Accept in args a list of tupple (sequence, name) and optionally the framesize (int)
    sequence_list = [(sequence, name), (sequence, name),...]
    """
    average_gc = [average_gc_by_frame(sequence, len_frame) for sequence, name_seq in sequence_list]
    name_seq = [name for sequence, name in sequence_list]

    # Create a histogram with for each frame in the sequence the average of GC
    fig = ff.create_distplot(average_gc, name_seq)
    plot(fig)

    # Create a histogram with for each average of GC in the sequence
    # the number of sequence with this average
    # <!> Work in progress, go around please and pay attention to your head when you going out <!>
    # <!> Experimental code : he look forward to trying <!>
    hist_data = []
    for data, name in zip(average_gc, name_seq):
        hist_dat, bin_pos = histogram(data, bins=100)
        """ Warning,
        the function histogram return 2 lists,
        a list with the number of seq in each bins
        and a list who delimit the position of each bins
        """
        x_start, x_stop, x_step = bin_pos[0], bin_pos[-1], bin_pos[-1]-bin_pos[-2]
        trace = go.Histogram(
            x=hist_dat,
            name=name_seq,
            xbins=dict(start=x_start, end=x_step, size=x_step),
            opacity=0.5)
        hist_data.append(trace)
    layout = go.Layout(
        title="Number of sequence by average GC",
        xaxis=dict(title="average GC"),
        yaxis=dict(title="number of sequence")
        bargap=0.2,
        bargroupgap=0.3
    )
    fig2 = go.Figure(hist_data, layout)
    plot(fig2)


#   ==================
#   the area of shame:
#   ==================


# ===== abandoned approach =====
# def bio_package_gc(sequence):
#   return GC(sequence)


# ===== abandoned approach =====
# def average_gc(sequence):
#     sequence = sequence.upper()
#     gc_count = 0
#     for nucl in sequence:
#         if nucl in ["G", "C", "S"]:
#             gc_count += 1
#     average = (gc_count / len(sequence))*100
#     return average


# def gc_deviating(sequence):
#     average_gc = average_gc_by_frame(sequence)
#     copy_average_gc = list(average_gc)
#     copy_average_gc.sort
#     seq_under_pourcentil = copy_average_gc[:int(len(average_gc)*0.025)]
#     seq_over_pourcentil = copy_average_gc[-int(len(average_gc)*0.025):]
#     seq_inside_pourcentil = [i for i in range]
#     return(seq_under_pourcentil,average_gc,seq_over_pourcentil)
