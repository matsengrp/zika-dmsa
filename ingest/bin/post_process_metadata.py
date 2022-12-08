#! /usr/bin/env python3

"""Reformat pandas DataTables for a pathogen build.

Expecting one argument, a NCBI Virus metadata.tsv file
"""
# ===== Dependencies
import argparse
import os
import sys

import numpy as np
import pandas as pd


def parse_args():
    # Main help command
    parser = argparse.ArgumentParser(
        description="Reformat a NCBI Virus metadata.tsv file for a pathogen build."
    )
    # Add first argument
    parser.add_argument(
        "--metadata", help="NCBI Virus metadata.tsv file.", required=True
    )
    parser.add_argument(
        "--outfile",
        help="Output file name, e.g. processed_metadata.tsv.",
        required=True,
    )

    return parser.parse_args()


# === Private methods


def _set_strain_name(record):
    """Check Isolate_s and Strain_s to find the strain name"""
    if record["strain"] != record["accession"]:
        strain_name = record["strain"]
    elif record["strain"] == record["accession"] and pd.notna(record["strain_s"]):
        strain_name = record["strain_s"]
    else:
        strain_name = record["strain"]

    return (
        strain_name.replace(" ", "_")
        .replace("-", "_")
        .replace(".", "_")
        .replace("(", "_")
        .replace(")", "_")
    )


def _set_url(record):
    """Set url column from accession"""
    return "https://www.ncbi.nlm.nih.gov/nuccore/" + str(record["accession"])


def _set_paper_url(record):
    """Set paper_url from publication"""
    if pd.notna(record["publications"]):
        paper_url = "https://www.ncbi.nlm.nih.gov/pubmed/" + str(record["publications"])
        return paper_url.split(",")[0]
    return ""


# === Main Method
def main():
    args = parse_args()
    df = pd.read_csv(args.metadata, sep="\t", header=0)

    # Mutate commands
    df["strain"] = df.apply(_set_strain_name, axis=1)
    df["url"] = df.apply(_set_url, axis=1)
    df["paper_url"] = df.apply(_set_paper_url, axis=1)
    df["authors"] = df["abbr_authors"]
    df["city"] = df["location"]

    # Format output
    METADATA_COLUMNS = [
        "strain",
        "accession",
        "genbank_accession_rev",
        "date",
        "updated",
        "region",
        "country",
        "division",
        "city",
        "authors",
        "url",
        "title",
        "paper_url",
    ]
    df.to_csv(args.outfile, sep="\t", index=False, columns=METADATA_COLUMNS)


if __name__ == "__main__":
    main()
