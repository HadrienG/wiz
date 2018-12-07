#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function


import re
import logging

import requests

from tqdm import tqdm
from requests_html import HTMLSession

from Bio import Entrez

Entrez.tool = "wiz"
Entrez.email = ""


class Assembly(object):
    """object representation of an assembly"""

    def __init__(self, entrez_dict):
        super(Assembly, self).__init__()
        self._all = entrez_dict
        self.organism = entrez_dict["Organism"]
        self.species_name = entrez_dict["SpeciesName"]

        self.taxid = entrez_dict["Taxid"]
        self.species_taxid = entrez_dict["SpeciesTaxid"]

        self.accession = entrez_dict["AssemblyAccession"]
        self.versioned_accession = entrez_dict["LastMajorReleaseAccession"]

        self.status = entrez_dict["AssemblyStatus"]
        self.name = entrez_dict["AssemblyName"]

        self.ftp_refseq = entrez_dict["FtpPath_RefSeq"]
        self.ftp_genbank = entrez_dict["FtpPath_GenBank"]


def query_assemblies(organism, output):
    """from a taxid or a organism name, download all refseq assemblies
    """
    logger = logging.getLogger(__name__)

    assemblies = []

    entrez_genomes = Entrez.read(Entrez.esearch(
        "assembly",
        term="%s[Organism]" % organism, retmax=10000))["IdList"]
    logger.info(f"Found {len(entrez_genomes)} organisms in ncbi assemblies")

    logger.info("Querying information about the assemblies. Be patient.")
    for id in tqdm(entrez_genomes):
        try:
            entrez_assembly = Entrez.read(
                Entrez.esummary(
                    db='assembly',
                    id=id))["DocumentSummarySet"]["DocumentSummary"][0]
        except KeyError as e:
            entrez_assembly = Entrez.read(
                Entrez.esummary(db='assembly', id=id))["DocumentSummarySet"]
            print(entrez_assembly.keys())
            raise
        else:
            assembly = Assembly(entrez_assembly)
            output_file = f"{output}/{assembly.accession}.fasta"
            download(assembly.ftp_refseq, output_file)
            assemblies.append(assembly)

    return assemblies


def download(url, output, chunk_size=1024):
    """download an url
    """
    logger = logging.getLogger(__name__)

    if url.startswith("ftp://"):
        url = url.replace("ftp://", "https://")
    if url:
        request = requests.get(url, stream=True)

        with open(output, "wb") as f:
            for chunk in request.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    f.flush()
