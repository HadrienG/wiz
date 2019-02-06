from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC


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
        error = "The sequence is void"
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
            average_gc = 0
            if pos+len_frame < len(sequence):
                average_gc = GC(sequence[pos:pos+len_frame])
            else:
                average_gc = GC(sequence[pos:])
            average_gc_by_windows.append(average_gc)
        return average_gc_by_windows


def gc_histogram(sequence, len_frame=5000):
    # coming soon
    pass
