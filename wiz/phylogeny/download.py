#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function


import re
import csv
import logging

import requests

from tqdm import tqdm

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


def query_assemblies(organism, output, quiet=False, representative=False):
    """from a taxid or a organism name, download all refseq assemblies
    """
    logger = logging.getLogger(__name__)

    assemblies = []

    if representative:
        genomes = Entrez.read(Entrez.esearch(
            "assembly",
            term=f"{organism}[Organism] AND \"representative genome\"[filter]",
            retmax=10000))["IdList"]
    else:
        genomes = Entrez.read(Entrez.esearch(
            "assembly",
            term=f"{organism}[Organism]",
            retmax=10000))["IdList"]
    logger.info(
        f"Found {len(genomes)} organisms in ncbi assemblies for {organism}")

    logger.info("Downloading the assemblies. Please be patient.")
    for id in tqdm(genomes, disable=quiet):
        try:
            entrez_assembly = Entrez.read(
                Entrez.esummary(
                    db="assembly",
                    id=id))["DocumentSummarySet"]["DocumentSummary"][0]
        except KeyError as e:
            entrez_assembly = Entrez.read(
                Entrez.esummary(db="assembly", id=id))["DocumentSummarySet"]
            print(entrez_assembly.keys())
            raise
        else:
            assembly = Assembly(entrez_assembly)
            output_file = f"{output}/{assembly.accession}.fasta"
            download(assembly.ftp_refseq, output_file)
            assemblies.append(assembly)

    return assemblies


def download(url, output_file, chunk_size=1024):
    """download an url
    """
    if url.startswith("ftp://"):  # requests doesnt support ftp
        url = url.replace("ftp://", "https://")
    if url:
        request = requests.get(url, stream=True)

        with open(output_file, "wb") as f:
            for chunk in request.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    f.flush()


def create_summary(assemblies, output_file):
    """create a summary csv file from a list of Assemblies
    """
    with open(output_file, "w", newline="") as csv_file:
        writer = csv.writer(csv_file, delimiter=",",
                            quotechar="|", quoting=csv.QUOTE_MINIMAL)
        header = ["Accession", "Organism", "Species", "Taxid", "Species Taxid"]
        writer.writerow(header)
        for assembly in assemblies:
            row = [assembly.accession,
                   assembly.organism,
                   assembly.species_name,
                   assembly.taxid,
                   assembly.species_taxid]
            writer.writerow(row)
