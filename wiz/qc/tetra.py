#!/usr/bin/env python
# -*- coding: utf-8


def tetranuc_count(sequence):
    tetra_dic = {}
    total_count = 0
    buffer = str(sequence[:3])
    for nucl in sequence[3:]:
        buffer += str(nucl)
        total_count +=1
        if buffer in tetra_dic.keys():
            tetra_dic[buffer] += 1
        else:
            tetra_dic[buffer] = 1
        buffer = str(buffer[1:])
    key_list = [key for key in tetra_dic.keys()]
    for key in key_list:
        tetra_dic[key] = (tetra_dic[key]/total_count)
    return tetra_dic


# * It's ok here
# TODO Make dataframe with the values returned
