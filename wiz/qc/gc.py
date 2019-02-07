from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from plotly.offline import plot
import plotly.figure_factory as ff


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


def gc_hist(sequence, len_frame=5000):
    data = average_gc_by_frame(sequence, len_frame)
    fig = ff.create_distplot(data, "bin", bin_size=100)
    plot(fig)

# def gc_deviating(sequence):
#     average_gc = average_gc_by_frame(sequence)
#     copy_average_gc = list(average_gc)
#     copy_average_gc.sort
#     seq_under_pourcentil = copy_average_gc[:int(len(average_gc)*0.025)]
#     seq_over_pourcentil = copy_average_gc[-int(len(average_gc)*0.025):]
#     seq_inside_pourcentil = [i for i in range]
#     return(seq_under_pourcentil,average_gc,seq_over_pourcentil)
