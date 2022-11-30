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
    description = "USAGE: python post_process_metadata.py --metadata 'data/metadata.tsv'"
  )
  # Add first argument
  parser.add_argument(
    "--metadata",
    help = "NCBI Virus metadata.tsv file.",
    required = True
  )
  parser.add_argument(
    "--outfile", default="processed_metadata.tsv",
    help = "Processed metdata file name. [default: processed_metadata.tsv] ",
    required = False
  )
  
  return parser.parse_args()

# === Private methods

def _set_strain_name(record):
  """Check Isolate_s and Strain_s to find the strain name"""
  strain_name=""
  if record['strain'] != record['accession']:
    strain_name=record['strain']
  elif record['strain'] == record['accession'] and pd.notna(record['strain_s']):
    strain_name=record['strain_s']
  else:
    strain_name=record['strain']

  return strain_name.replace(' ', '_').replace('-',"_").replace(".","_").replace("(","_").replace(")","_")

def _set_url(record):
  """Set url column from accession"""
  return "https://www.ncbi.nlm.nih.gov/nuccore/" + record['accession']

def _set_paper_url(record):
  """Set paper_url from publication"""
  if pd.notna(record['publications']):
    return "https://www.ncbi.nlm.nih.gov/pubmed/" + record['publications']
  return ""

def _set_dengue_serotype(record):
  """Set dengue serotype from viruslineage_ids"""
  if "11053" in record['viruslineage_ids']:
    return "denv1"
  elif "11060" in record['viruslineage_ids']:
    return "denv2"
  elif "11069" in record['viruslineage_ids']:
    return "denv3"
  elif "11070" in record['viruslineage_ids']:
    return "denv4"
  else:
    return ""

# === Main Method
def main():
  args = parse_args()
  df=pd.read_csv(args.metadata, sep='\t', header=0)

  # Mutate commands
  df['strain'] = df.apply(_set_strain_name, axis=1)
  df['url'] = df.apply(_set_url, axis=1)
  df['paper_url'] = df.apply(_set_paper_url, axis=1)
  df['serotype'] = df.apply(_set_dengue_serotype, axis=1)
  df['authors'] = df['abbr_authors']
  df['city'] = df['location']

  # Format output
  METADATA_COLUMNS = ["strain", "accession", "genbank_accession_rev","serotype", "date", "updated", "region", "country", "division", "city", "authors", "url", "title", "paper_url"]
  df = df.reindex(METADATA_COLUMNS, axis=1)
  df.to_csv(args.outfile, sep='\t', index=False)

if __name__ == '__main__':
  main()