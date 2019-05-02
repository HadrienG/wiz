#!/usr/bin/env python
# -*- coding: utf-8

"""
File that gathers the necessary functions for taxonomics operations.
"""

import json
import os
import logging

from taxadb.accessionid import AccessionID
from taxadb.names import SciName
from taxadb.taxid import TaxID

from wiz.qc import tools


logger = logging.getLogger(__name__)


def taxonomy(g_id, g_seq, args):
    """
    main function of module taxonomy
    """
    logger.info(f"Sketching {g_id}")
    extract_contig(g_id, g_seq, args)
    tools.finch_sketch(g_id, args.output)
    logger.info(f"Search for sequences similar to {g_id}")
    tools.finch_dist(g_id, args.finchdb, args.output)
    logger.info(f"Identification of the taxonomy of the results")
    identification = extract_finch_data(g_id, args)
    tax = taxid_search(identification, args)
    return tax


def taxid_search(dist_results, args):
    """
    Function who's looking for tax id based on a UID or a scientitific name.
    """
    # taxadb link init
    names = SciName(dbtype='sqlite', dbname=args.taxadb)
    accession = AccessionID(dbtype='sqlite', dbname=args.taxadb)
    taxid = TaxID(dbtype='sqlite', dbname=args.taxadb)
    tax = []
    for result in dist_results:
        reference, jaccard = result
        if reference != "Nothing found":
            tax_id, lineage = 0, []
            tax_id = accession.taxid(reference)
            if not isinstance(tax_id, int):  # type(tax_id) != int:
                reference = reference.split(".")[0]
                tax_id = accession.taxid(reference)
            if not isinstance(tax_id, int):  # type(tax_id) != int:
                tax_id = names.taxid(reference)
            if isinstance(tax_id, int):  # type(tax_id) == int:
                lineage = taxid.lineage_name(tax_id, reverse=False, ranks=True)
                logger.debug(f" Taxid :   found   : {reference}.")
                tax.append((lineage, jaccard))
            else:
                logger.debug(f" Taxid : not found : {reference}.")
                tax.append((
                    {"species": reference, "no rank": "Not found"}, jaccard))
        else:
            logger.info(f" Taxid : not found : {reference}.")
            tax.append((
                {"species": reference, "no rank": "Not found"}, jaccard))
    return tax


def extract_contig(contig_id, contig, args):
    """
    Function who's write contig in a file in fasta format.
    """
    with open(f"{args.output}/finch/{contig_id}.fna", "w") as file:
        contig = str(contig)
        file.write(f"> {contig_id}\n")
        for i in range(0, len(contig), 80):
            file.write(contig[i:i+80]+"\n")
        file.write("\n\n")


def extract_finch_data(filename, args):
    """
    Function who's parse data from the finch output file.
    """
    finch_out = []
    finch_file = open(f"{args.output}/finch/{filename}.json", "r").readlines()
    for line in finch_file:
        results = json.loads(line)
        if results:
            for result in results:
                finch_out.append((result["reference"], result["jaccard"]))
    logger.info(f" {len(finch_out)} results found for {filename}")
    if not finch_out:  # len(finch_out) == 0:
        finch_out.append(("Nothing found", 1))
    os.remove(f"{os.path.abspath(args.output)}/finch/{filename}.fna")
    os.remove(f"{os.path.abspath(args.output)}/finch/{filename}.sk")
    os.remove(f"{os.path.abspath(args.output)}/finch/{filename}.json")
    # return a list of tuple with each result of distance calculation
    # finch_out=[(ref, dist), (ref, dist), ...]
    return finch_out
