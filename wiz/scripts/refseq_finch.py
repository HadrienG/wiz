#!/usr/bin/env python3
# -*- coding: utf-8

import os
import subprocess


folder_count, finch = 0, 1
tot_folder = len(os.listdir("/db/refseq/bacteria/folder/"))
for i in range(0, tot_folder-1):
    if not os.path.isfile(f"/db/finch_db/bacteria_N_k21.part{finch}.sk"):
        print(f"sketching folder_{i}")
        s_path = f"/db/refseq/bacteria/folder/folder_{i}"
        seqs_name = os.listdir(s_path)
        seqs = [f"{s_path}/{name}" for name in seqs_name]
        arg = ["finch", "sketch", "-N", "-o", f"/db/finch_db/bacteria_N_k21.part{finch}.sk"]
        arg += seqs
        subprocess.run(arg)
    else:
        print(f"folder_{i} ok!")
    folder_count += 1
    finch += 1
