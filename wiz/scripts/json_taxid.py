#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Use taxadb to find the taxonomic id of sequence in a fasta file on the basis of their name and save it all in a json file, if the program does not find the identifier the user can enter it manually. If the user does not validate a number and press enter the program will consider that the identifier was not found.
"""
from Bio import SeqIO
import argparse
import os
from os import path
import json

parser = argparse.ArgumentParser(
    prog="refseq dl",
    usage="refseq [options] ...",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument(
    "-f",
    "--file",
    type=str,
    metavar="",
    help="input file in fasta format",
    required=True,
    nargs='+'
)
parser.add_argument(
    "-o",
    "--out",
    type=str,
    metavar="",
    help="output file",
    required=True,
    #nargs='+'
)
parser.add_argument(
    "-t",
    "--taxadb",
    type=str,
    metavar="",
    help="sqlite taxadb file",
    required=True,
    #nargs='+'
)


args = parser.parse_args()

count_ok, count_fail, count = 0, 0, 0
from taxadb.names import SciName
print("loading db")
names = SciName(dbtype='sqlite', dbname=args.taxadb)
total_seq_count = 0
million_seq = 0
models = {}
global_Seq = {}
print("checking backup")
if os.path.isfile(f"{args.out}_temp.json"):
    print(f"{args.out}_temp.json exist, it used to recovered id")
    with open(f"{args.out}_temp.json", "r") as f:
        global_Seq = json.load(f)
else:
    print("no backup detected")

print(args.out)
for file in args.file:
    fname = path.basename(file)
    fasta_file = SeqIO.parse(file, "fasta")
    print(f"\033[92mspliting {file} \033[0m")
    model = fname.split(".")[0]
    seq = {}
    for record in fasta_file:
        desc = record.description
        if "{" in desc:
            name = str(desc[desc.index("{")+1:desc.index("}")])
        else:
            bracket_count = 0
            for i in desc:
                if i == "[":
                    bracket_count += 1
            if bracket_count == 1:
                name = str(desc[desc.index("[")+1:desc.index("]")])
            else:
                desc_temp = list(desc)
                for i in range(0, bracket_count-1):
                    desc_temp.pop(desc_temp.index("["))
                    desc_temp.pop(desc_temp.index("]"))
                desc_temp2 = ""
                for i in desc_temp:
                    desc_temp2 += i
                name = str(desc_temp2[desc_temp2.index("[")+1:desc_temp2.index("]")])
            if "," in name:
                name = str(name[:name.index(",")])
        blk = ' '*8
        print(f"search :{blk}{name}", end="")
        total_seq_count += 1
        if name in global_Seq.keys():
            seq[desc] = global_Seq[name]
            tax = global_Seq[name] # line for print info
            print(f"\rsearch : \033[32m{tax:^8}\033[0m {name}")
        else:
            tax = names.taxid(name)
            # if not path.isdir(f"{args.out}"):
            #     os.makedirs(f"{args.out}/ok")
            #     os.makedirs(f"{args.out}/fail")
            if type(tax) == int:
                print(f"\rsearch : \033[34m{tax:^8}\033[0m {name}")
                # print('\033[34m', tax, '\033[0m')
                seq[desc] = tax

            else:
                fail = "FAIL"
                print(f"\rsearch : \033[38;5;166m{fail:^8}\033[0m {name}")
                # print('\033[91m', "FAIL", '\033[0m')
                rep = ""
                print(f"DESC : {desc}")
                os.system(f"echo '{name}' | xclip -selection clipboard")
                with open(f"{args.out}_temp.json", "w") as f:
                    backup = json.dumps(global_Seq)
                    f.writelines(backup)
                while rep != "y":
                    tax = input("Taxid : ")
                    if tax == "":
                        print("tax id = \033[91mNo found\33[0m? (y/N)")
                        taxid = "No found"
                        rep = input()
                    else:
                        try:
                            tax = int(tax)
                            print(f"tax id = \033[94m{tax}\33[0m (y/N)")
                            rep = input()
                        except:
                            print(f"tax id = \033[91mTEXT\33[0m")
                            print("\033[91minteger value only\33[0m")
                            rep = ""
                seq[desc] = tax
            global_Seq[name] = tax
    with open(f"{args.out}_temp.json", "w") as f:
        backup = json.dumps(global_Seq)
        f.writelines(backup)
    models[model] = seq
    if len(models) % 100 == 0:
        to_write = json.dumps(models)
        with open(f"{args.out}.json", "w") as f:
            f.writelines(to_write)
    if total_seq_count >= 1000000:
        million_seq +=1
        total_seq_count -= 1000000
print("Final saving")
to_write = json.dumps(models)
with open(f"{args.out}.json", "w") as f:
    f.writelines(to_write)
print(f"models : {len(models)}")
print(f"sequences tot : {million_seq}{total_seq_count}")
print(f"dif seq: {len(global_Seq)} ")
print("done")
