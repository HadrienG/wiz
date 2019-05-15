# graphs.py

def distplot_gc(contigs):
    # TODO fix bug here, if a contig have a only one value for averages_gc
    this next line crash_ 
    line 188 
    fig = distplot(gc_averages, seq_names, show_hist=False)


def get_nodes()
line 458 to 478 !
# TODO consider the case where Jaccard's distance is 0. 1 divide by jaccard will cause a bug 

error with gc bounds in scatter_GC:
    the bounds must be for global GC_content of a bin and not for each contigs inside the bin
    