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
    logger.debug(f" Coding density : {((coding_region-overlap_region)/len_seq)*100}")
    return ((coding_region-overlap_region)/len_seq)*100


def taxonomy(g_id, g_seq, args):
    logger.info(f"Sketching {g_id}")
    extract_contig(g_id, g_seq, args)
    tools.finch_sketch(g_id, args.output)
    logger.info(f"Search for sequences similar to {g_id}")
    tools.finch_dist(g_id, args.finchdb, args.output)
    logger.info(f"Identification of the taxonomy of the results")
    identification = extract_finch_data(g_id, args)
    tax = taxid(identification, args)
    return tax


def taxid(dist_results, args):
    # taxadb link init
    names = SciName(dbtype='sqlite', dbname=args.taxadb)
    accession = AccessionID(dbtype='sqlite', dbname=args.taxadb)
    taxid = TaxID(dbtype='sqlite', dbname=args.taxadb)
    tax = []
    for result in dist_results:
        ref, jaccard = result
        if ref != "Nothing found":
            tax_id, lineage = 0, []
            tax_id = accession.taxid(ref)
            if type(tax_id) != int:
                r = ref.split(".")[0]
                tax_id = accession.taxid(r)
            if type(tax_id) != int:
                tax_id = names.taxid(ref)
            if type(tax_id) == int:
                lineage = taxid.lineage_name(tax_id, reverse=False, ranks=True)
                logger.debug(f" Taxid :   found   : {ref}.")
                tax.append((lineage, jaccard))
            else:
                logger.debug(f" Taxid : not found : {ref}.")
                tax.append(({"species": ref, "no rank": "Not found"}, jaccard))
        else:
            logger.info(f" Taxid : not found : {ref}.")
            tax.append(({"species": ref, "no rank": "Not found"}, jaccard))
    return tax


def extract_contig(seq_id, seq, args):
    with open(f"{args.output}/finch/{seq_id}.fna", "w") as fw:
        seq = str(seq)
        fw.write(f"> {seq_id}\n")
        for i in range(0, len(seq), 80):
            fw.write(seq[i:i+80]+"\n")
        fw.write("\n\n")


def extract_finch_data(filename, args):
    finch_out = []
    finch_file = open(f"{args.output}/finch/{filename}.json", "r").readlines()
    for line in finch_file:
        result = json.loads(line)
        if result != 0:
            for r in result:
                finch_out.append((r["reference"], r["jaccard"]))
    logger.info(f" {len(finch_out)} results found for {filename}")
    if len(finch_out) == 0:
        finch_out.append(("Nothing found", 0.5))
    os.remove(f"{os.path.abspath(args.output)}/finch/{filename}.fna")
    os.remove(f"{os.path.abspath(args.output)}/finch/{filename}.sk")
    os.remove(f"{os.path.abspath(args.output)}/finch/{filename}.json")
    # return a list of tuple with each result of distance calculation
    # finch_out=[(ref, dist), (ref, dist), ...]
    return(finch_out)
