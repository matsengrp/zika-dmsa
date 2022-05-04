#! /usr/bin/env bash

set -v

# TODO: Generalizable flags: family-species, fromyear, minlength, metadata

curl "https://www.viprbrc.org/brc/api/sequence?datatype=genome&family=flavi&species=Zika%20virus&fromyear=2013&minlength=5000&metadata=genbank,strainname,segment,date,host,country,genotype,species&output=fasta" | \
  tr '-' '_' | \
  tr ' ' '_' | \
  sed 's:N/A:NA:g' > \
  GenomicFastaResults.fasta

echo "total sequences: " `grep -c ">" GenomicFastaResults.fasta`
