#!/usr/bin/env python
# -*- coding: utf-8


def tetranuc_count(sequence):
    tetra_dic = {}
    buffer = str(sequence[:3])
    for nucl in sequence[3:]:
        buffer += str(nucl)
        if buffer in tetra_dic.keys():
            tetra_dic[buffer] += 1
        else:
            tetra_dic[buffer] = 1
        buffer = str(buffer[1:])
    return tetra_dic