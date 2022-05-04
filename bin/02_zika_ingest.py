#! /usr/bin/env python3

# ==== Packages
import os, re, time, datetime, csv, sys
from Bio import SeqIO
from typing import NamedTuple

# === Input variables
zika_seqs = "../example_data/small.fasta"

# From fauna
strain_fix_fname =  "zika_strain_name_fix.tsv"
location_fix_fname = "zika_location_fix.tsv"
date_fix_fname = "zika_date_fix.tsv"

virus_fasta_fields = {1:'strain', 3:'collection_date', 4: 'host', 5:'country'}
sequence_fasta_fields = {0:'accession', 1:'strain'}

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

# Originally function: replace_strain_name, also at L112 in upload.py
def fixes_str(original_str:str, fixes_dict:dict[str,str]={}):
    '''
    return the new strain name/location/date that will replace the original string
    '''
    if original_str in fixes_dict:
        return fixes[original_str] 
    else:
        return original_str

# vdl/zika_uploads.py
def fixes_strain_name(name, fixes_tsv:str="") -> (str,str): # Since we can't decide if we want strain or name
    fixes_dict = {}
    if(len(fixes_tsv)>0):
        fixes_dict = define_fixes_dict(fixes_tsv)
    
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
    return name, original_name

# vdl/parse.py Load data
def parse_fasta_file(fasta, virus_fasta_fields, sequence_fasta_fields, **kwargs):
    '''
    Parse FASTA file with default header formatting
    :return: list of documents(dictionaries of attributes) to upload
    '''
    header_fixes = False
    if (kwargs["fasta_header_fix"]):
        header_fixes = {}
        try:
            with open(kwargs["fasta_header_fix"], 'rU') as fh:
                for line in fh:
                    if not line.startswith('#'):
                        k, v = line.strip().split("\t")
                        header_fixes[k] = v                
        except IOError:
            raise Exception(kwargs["fasta_header_fix"], "not found")
    viruses = []
    sequences = []
    try:
        handle = open(fasta, 'r')
    except IOError:
        raise Exception(fasta, "not found")
    else:
        for record in SeqIO.parse(handle, "fasta"):
            if header_fixes:
                try:
                    record.description = header_fixes[record.description]
                except KeyError:
                    raise Exception(record.description, "not in header fix file. Fatal.")
            content = list(map(lambda x: x.strip(), record.description.replace(">", "").split('|')))
            v = {key: content[ii] if ii < len(content) else "" for ii, key in virus_fasta_fields.items()}
            s = {key: content[ii] if ii < len(content) else "" for ii, key in sequence_fasta_fields.items()}
            s['sequence'] = str(record.seq).lower()
            #v = self.add_virus_fields(v, **kwargs)
            #s = self.add_sequence_fields(s, **kwargs)
            sequences.append(s)
            viruses.append(v)
        handle.close()
    return (viruses, sequences)

# === Only fix casing on the Host?
def fix_casing(self, document): # JC
    for field in ['host']:
        if field in document and document[field] is not None:
            document[field] = self.camelcase_to_snakecase(document[field])

### === Main Method
# vdl/upload.py L58
# self.connect(**kwargs)
# print("Uploading Viruses to VDB")
# viruses, sequences = self.parse(**kwargs) #<= expand vdl/parse.py
zika_fasta = parse_fasta_file(zika_seqs, virus_fasta_fields, sequence_fasta_fields, fasta_header_fix = False)
viruses, sequences = parse_fasta_file(zika_seqs, virus_fasta_fields, sequence_fasta_fields, fasta_header_fix = False)
# print('Formatting documents for upload')
# self.format_viruses(viruses, **kwargs)
# self.format_sequences(sequences, **kwargs)

# print("")
# print("Filtering out viruses")
# viruses = self.filter(viruses, 'strain', **kwargs)
# print("Filtering out sequences")
# sequences = self.filter(sequences, 'accession', **kwargs)
# print("")

# === Test combined dictionary methods
fix_name_dict = define_fixes_dict(strain_fix_fname)  # tsv file in Input
fix_location_dict = define_fixes_dict(location_fix_fname)
fix_date_dict = define_fixes_dict(date_fix_fname)

# === Process a smaller data file
zika_fasta = "../example_data/small.fasta"

zika_seqs = parse_fasta_file(zika_seqs, virus_fasta_fields, sequence_fasta_fields, fasta_header_fix = False)

print(zika_seqs[0])
print("\n\nstrain: ",zika_seqs[0][0]['strain'])
print("fix_strain_name output:", fixes_strain_name(zika_seqs[0][0]['strain']))