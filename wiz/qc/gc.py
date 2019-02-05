from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC


def bio_package_gc(sequence):
    return GC(sequence)


def average_gc(sequence):
    sequence = sequence.upper()
    gc_count = 0
    for nucl in sequence:
        if nucl in ["G", "C", "S"]:
            gc_count += 1
    average = (gc_count / len(sequence))*100
    return average


def average_gc_by_frame(sequence, len_frame):
    if frame > len(sequence):
        raise ValueError("The size of the frame is superior of the sequence length.")
    if frame <= 0:
        raise ValueError("The size of the frame is negative or null.")
    if len(sequence) == 0:
        raise ValueError("The sequence is void")
    average_gc_by_windows = []
    for pos in range(0, len(sequence), len_frame):
        average_gc = 0
        if pos+frame < len(sequence):
            average_gc = GC(sequence[pos:pos+frame])
        else:
            average_gc = GC(sequence[pos:])
        average_gc_by_windows.append(average_gc)
    return average_gc_by_windows
