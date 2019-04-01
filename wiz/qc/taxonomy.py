#!/usr/bin/env python
# -*- coding: utf-8

from wiz.annotate.tools import prodigal
from wiz.qc import tools
from taxadb.accessionid import AccessionID
from taxadb.names import SciName
from taxadb.taxid import TaxID
import json
import os
import logging



logger = logging.getLogger(__name__)


def coding_density(genes_pos, len_seq):
    coding_region = 0
    overlap_region = 0
    for genes in genes_pos:
        coding_region += (genes[1]-genes[0])
    for genes in genes_pos[:-1]: 
        n_genes = genes_pos[genes_pos.index(genes)+1]
        overlap_region += max(0, min(genes[1], n_genes[1])-max(genes[0], n_genes[0]))
    # print(coding_region,overlap_region)
    # print(((coding_region-overlap_region)/len_seq)*100)
    return ((coding_region-overlap_region)/len_seq)*100


# ok to remove this function
def coding_density2(genes_pos, len_seq):
    # print(genes_pos)
    coding_region = [(genes_pos[0])]
    for genes in genes_pos[1:]:
        r_start, r_stop = coding_region[-1]
        g_start, g_stop = genes
        if g_start > r_stop+1:
            coding_region.append(genes)
        else:
            if g_stop > r_stop:
                coding_region[-1] = (r_start, g_stop)
    tot_coding = 0
    for region in coding_region:
        tot_coding += (region[1]-region[0])
    # print(tot_coding)
    tot_coding = (tot_coding/len_seq)*100
    # print(tot_coding)
    return tot_coding
# =========================


def taxonomy(dist_results, args, bins):
    # taxadb link init
    names = SciName(dbtype='sqlite', dbname=args.taxadb)
    accession = AccessionID(dbtype='sqlite', dbname=args.taxadb)
    taxid = TaxID(dbtype='sqlite', dbname=args.taxadb)
    for key in dist_results:
        tax_id = 0
        lineage = []
        results = dist_results[key]
        r_temp = []
        for result in results:
            ref, jaccard = result
            if isinstance(accession.taxid(ref), int):
                tax_id = accession.taxid(ref)
            else:
                taxid = names.taxid(ref)
            if type(tax_id) == int:
                lineage = taxid.lineage_name(tax_id)
                r_temp.append(ref, lineage, jaccard)
            else:
                raise ValueError(f"Unfound taxid for {ref}. Please check the databases !")
        dist_results[key] = r_temp
    for genome_bin in bins:
        if genome_bin.id in dist_results.keys():
            genome_bin.taxonomy = dist_results[genome_bin.id]


def extract_contig(seq_id, seq, args):
    with open(f"{args.output}/finch/{seq_id}.fna", "w") as fw:
        seq = str(seq)
        fw.write(f"> {seq_id}\n")
        for i in range(0, len(seq), 80):
            fw.write(seq[i:i+80]+"\n")
        fw.write("\n\n")


def sketching_contig(contigs, filename, args):
    filename = os.path.basename(filename)
    tools.finch_sketch(filename, contigs, args.output)
    #for contig in contigs:
    #    os.remove(f"{args.output}/finch/{contig}.fna")


def dist_contigs(filename, args):
    tools.finch_dist(filename, args.finchdb, args.output)
    finch_out = {}
    finch_file = open(f"{args.output}/finch/{filename}.json", "r").readlines()
    for line in finch_file:
        result = json.loads(line)
        if result != 0:
            for r in result:
                if r["query"] in finch_out.keys():
                    finch_out[r["query"]].append((r["reference"], r["jaccard"]))
                else:
                    finch_out[r["query"]] = [(r["reference"], r["jaccard"])]
    logger.info(f" {len(finch_out)} results found")
    os.remove(f"{os.path.abspath(args.output)}/finch/{filename}.sk")
    os.remove(f"{os.path.abspath(args.output)}/finch/{filename}.json")
    # return dictionnary with list of tuple with each result of distance calculation
    # finch_out[contig]=[(ref, dist), (ref, dist), ...]
    return(finch_out)
