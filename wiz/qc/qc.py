#!/usr/bin/env python
# -*- coding: utf-8

import gc


def run(args):
    """
    main function for wiz quality check (qc)
    """
    # speculative code, I don't think completely how to use the args.xxx function
    # args.bins = [(sequence,name),...]
    args.window = 5000
    args.hist_average_per_window = True
    args.hist_nbr_seq_per_average = True
    gc_per_bin = []
    for bin in args.bins:
        seq, name = bin
        average_gc = gc.average_gc(seq, args.window, truncate=False)
        gc_per_bin.append((average_gc, name))

    if args.hist_average_per_window:
        gc.scatter_gc(gc_per_bin, args.window)

    if args.hist_nbr_seq_per_average:
        gc.distplot_gc(gc_per_bin)
        gc.histogram_gc(gc_per_bin)

    print("implemented!")
