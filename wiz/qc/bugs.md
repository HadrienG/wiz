# graphs.py

def distplot_gc(contigs):
    _TODO fix bug here, if a contig have a only one value for averages_gc
    this next line crash_
    fig = distplot(gc_averages, seq_names, show_hist=False)
