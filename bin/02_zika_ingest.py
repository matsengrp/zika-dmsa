#! /usr/bin/env python3

# ==== Packages
import os, re, time, csv, sys
from Bio import SeqIO
from typing import NamedTuple
from datetime import datetime

# === Input variables
zika_fasta = "../example_data/small.fasta"

# From fauna
strain_fix_fname =  "zika_strain_name_fix.tsv"
location_fix_fname = "zika_location_fix.tsv"
date_fix_fname = "zika_date_fix.tsv"

#virus_fasta_fields = {1:'strain', 3:'collection_date', 4: 'host', 5:'country'}
#sequence_fasta_fields = {0:'accession', 1:'strain'}
# Seems duplicative, replace with:
header_fasta_fields = {0:'accession', 1:'strain', 3:'collection_date', 4: 'host', 5:'country'}
# If we're ignoring anything past 5, then don't pull from vipr

virus="zika" # upload
expected_date_formats = {'%Y_%m_%d','%Y (Month and day unknown)', '%Y-%m (Day unknown)', 
                         '%Y %b %d','%Y-XX-XX', '%Y-%m-%d', '%Y-%m-%dT%H:%M:%SZ'}

# ==== Functions
# vdl/uploads.py
# Originally 3 functions: define_strain_fixes, define_location_fixes, define_date_fixes
def define_fixes_dict(fname:str) -> dict[str,str]:
    '''
    Open strain/location/date fixing tsv files and define corresponding dictionaries
    '''
    reader = csv.DictReader(filter(lambda row: row[0]!='#', open(fname)), delimiter='\t')
    fixes_dict = {}
    for line in reader:
        fixes_dict[line['label'].encode().decode('unicode-escape')] = line['fix']
    return fixes_dict

def fixes_str(original_str:str, fixes_dict:dict[str,str]={}):
    '''
    return the new strain name that will replace the original string
    cannot be applied to location and date since key is based on strain name
    '''
    # labmda x: fixes[original_str] if original_str in fixes dict else original_str
    if original_str in fixes_dict:
        return fixes[original_str] 
    else:
        return original_str

# vdl/zika_uploads.py
def fixes_strain_name(name, fixes_dict={}, fixes_tsv:str="") -> (str,str): # Since we can't decide if we want strain or name
    if(len(fixes_dict)<1 and len(fixes_tsv)>0):
        fixes_dict = define_fixes_dict(fixes_tsv)
    
    # https://stackoverflow.com/questions/2484156/is-str-replace-replace-ad-nauseam-a-standard-idiom-in-python
    original_name = name
    name = fixes_str(original_name, fixes_dict) 
    name = name.replace('Zika_virus', '').replace('Zikavirus', '').replace('Zika virus', '').replace('Zika', '').replace('ZIKV', '')
    name = name.replace('Human', '').replace('human', '').replace('H.sapiens_wt', '').replace('H.sapiens_tc', '').replace('Hsapiens_tc', '').replace('H.sapiens-tc', '').replace('Homo_sapiens', '').replace('Homo sapiens', '').replace('Hsapiens', '').replace('H.sapiens', '')
    name = name.replace('/Hu/', '')
    name = name.replace('_Asian', '').replace('_Asia', '').replace('_asian', '').replace('_asia', '')
    name = name.replace('_URI', '').replace('_SER', '').replace('_PLA', '').replace('_MOS', '').replace('_SAL', '')
    name = name.replace('Aaegypti_wt', 'Aedes_aegypti').replace('Aedessp', 'Aedes_sp')
    name = name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('//', '/').replace('__', '_').replace('.', '').replace(',', '')
    name = re.sub('^[\/\_\-]', '', name)
    try: # ID must start with letter
        name = 'V' + str(int(name))
    except:
        pass
    name = fixes_str(name, fixes_dict)
    return name

def ncov_ingest_format_date(date_string: str, expected_formats: set) -> str:
    """
    Format *date_string* to ISO 8601 date (YYYY-MM-DD).
    If *date_string* does not match *expected_formats*, return *date_string*.
    >>> expected_formats = {'%Y-%m-%d', '%Y-%m-%dT%H:%M:%SZ'}
    >>> format_date("2020", expected_formats)
    '2020'
    >>> format_date("2020-01", expected_formats)
    '2020-01'
    >>> format_date("2020-1-15", expected_formats)
    '2020-01-15'
    >>> format_date("2020-1-1", expected_formats)
    '2020-01-01'
    >>> format_date("2020-01-15", expected_formats)
    '2020-01-15'
    >>> format_date("2020-01-15T00:00:00Z", expected_formats)
    '2020-01-15'
    >>> format_date("2020-XX-XX", expected_formats)
    '2020'
    >>> format_date("2020 (Month and day unknown)", expected_formats)
    '2020'
    >>> format_date("2020-06 (Day unknown)", expected_formats)
    '2020-06'
    """
    for date_format in expected_formats:
        try:
            return datetime.strptime(date_string, date_format).strftime('%Y-%m-%d')
        except ValueError:
            continue
    
    return date_string

### === Main Method
# vdl/upload.py L58
# self.connect(**kwargs)
# print("Uploading Viruses to VDB")
# viruses, sequences = self.parse(**kwargs) #<= expand vdl/parse.py
# zika_fasta = parse_fasta_file(zika_seqs, virus_fasta_fields, sequence_fasta_fields, fasta_header_fix = False)
# viruses, sequences = parse_fasta_file(zika_seqs, virus_fasta_fields, sequence_fasta_fields, fasta_header_fix = False)
# print('Formatting documents for upload')
# self.format_viruses(viruses, **kwargs)
# self.format_sequences(sequences, **kwargs)

# print("")
# print("Filtering out viruses")
# viruses = self.filter(viruses, 'strain', **kwargs)
# print("Filtering out sequences")
# sequences = self.filter(sequences, 'accession', **kwargs)
# print("")

fname=zika_fasta

# Early exit if file not found
try:
    fhandle = open(fname, 'r')
except IOError:
    raise Exception(fname, "not found")

fix_name_dict = define_fixes_dict(strain_fix_fname)         # if defined 
fix_location_dict = define_fixes_dict(location_fix_fname)
fix_date_dict = define_fixes_dict(date_fix_fname)

try:
    shandle = open("sequences.fasta", 'w')
    mhandle = open("metadata.tsv", 'w')
    shandle.close()
    mhandle.close()
    shandle = open("sequences.fasta", 'a')
    mhandle = open("metadata.tsv", 'a')
    mhandle.write("\t".join(("strain", "virus", "accession", "date", "region", "country", "division", "city", "db", "segment", "authors", "url", "title", "journal", "paper_url"))+"\n")
except IOError:
    raise Exception('Cannot write to sequences.fasta and/or metadata.tsv')

for record in SeqIO.parse(fhandle, "fasta"):
    #print(record.id) # Header, Breaks at spaces!
    print(record.description) # Whole header
    content = list(
        map(lambda x: x.strip(), 
            record.description
            .replace(" ", "_") # Deal with spaces
            .split('|'))
    )
    print(content)
    metadata = {key: content[ii] if ii < len(content) else "" for ii, key in header_fasta_fields.items()}
    print("metadata=", metadata)
    metadata["strain"] = fixes_strain_name(metadata["strain"], fixes_dict=fix_name_dict)
    # metadata["collection_date"] = fixes_str(metadata["strain"], fix_date_dict) # based on fixed strain name?
    # metadata["location"] = fixes_str(metadata["strain"], fix_location_dict)  # find where location is defined
    print("metadata[strain]=", metadata["strain"])
    
    # Hmm, was an obj method (checking for "date", "collection date", "submission date", seems too specialized...)
    # If you want to check for all dates, then check all fields for a XXXX-XX-XX or similar date format...
    
    metadata["collection_date"]=ncov_ingest_format_date(metadata["collection_date"], expected_date_formats)
    print("metadata[collection_date]=", metadata["collection_date"])
    shandle.write(">" + metadata["strain"] + "\n")
    shandle.write(str(record.seq).lower() + "\n") # todo: split in lines of 80 char
    mhandle.write("\t".join((metadata["strain"], "zika", metadata["accession"], metadata["collection_date"]))+"\n")

shandle.close()
mhandle.close()
fhandle.close()
