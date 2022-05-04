### README

Summary of scripts and files. The numbers are added for steps of the process but will be removed in the final version.

```
bin
  |_ 01_pull_vipr.sh     #<= todo: mv to python and read directly into SeqIO
  |_ 02_zika_ingest.py   #<= should have seq obj, filters or fixes applied to header/seqs, subsets applied to headers/seqs
  |
  | # Fixes dictionaries for 02_zika_ingest.py # todo: combine into one look up table? On the fence on this...
  |_ zika_date_fix.tsv
  |_ zika_location_fix.tsv
  |_ zika_strain_name_fix.tsv
  |
  | # historical diagnostic script?  
  |_ check-countries-have-colors.sh
```
