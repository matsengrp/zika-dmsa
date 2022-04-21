#! /usr/bin/env python3

# ==== Packages
import os, re, time, datetime, csv, sys
from Bio import SeqIO

# === INput variables

zika_seqs = "../example_data/small.fasta"

# From fauna, think how to move over
strain_fix_fname =  "zika_strain_name_fix.tsv"
location_fix_fname = "zika_location_fix.tsv"
date_fix_fname = "zika_date_fix.tsv"
virus_fasta_fields = {1:'strain', 3:'collection_date', 4: 'host', 5:'country'}
sequence_fasta_fields = {0:'accession', 1:'strain'}

# ==== Functions
# vdl/uploads.py
def replace_strain_name(original_name, fixes={}):
    '''
    return the new strain name that will replace the original
    '''
    if original_name in fixes:
        return fixes[original_name] 
    else:
        return original_name 

# vdl/uploads.py
def define_strain_fixes(fname):  
    '''
    Open strain name fixing files and define corresponding dictionaries
    '''
    reader = csv.DictReader(filter(lambda row: row[0]!='#', open(fname)), delimiter='\t')
    fix_whole_name = {}
    for line in reader:
        fix_whole_name[line['label'].encode().decode('unicode-escape')] = line['fix']
    return fix_whole_name         

fix_whole_name = define_strain_fixes(strain_fix_fname)  # tsv file in Input
#type(fix_whole_name)
#print(list(fix_whole_name)[:10])

def fix_name(name): # polymorphism, overwrite fix_names in upload.py
    original_name = name
    name = replace_strain_name(original_name, fix_whole_name) 
    name = name.replace('Zika_virus', '').replace('Zikavirus', '').replace('Zika virus', '').replace('Zika', '').replace('ZIKV', '')
    name = name.replace('Human', '').replace('human', '').replace('H.sapiens_wt', '').replace('H.sapiens_tc', '').replace('Hsapiens_tc', '').replace('H.sapiens-tc', '').replace('Homo_sapiens', '').replace('Homo sapiens', '').replace('Hsapiens', '').replace('H.sapiens', '')
    name = name.replace('/Hu/', '')
    name = name.replace('_Asian', '').replace('_Asia', '').replace('_asian', '').replace('_asia', '')
    name = name.replace('_URI', '').replace('_SER', '').replace('_PLA', '').replace('_MOS', '').replace('_SAL', '')
    name = name.replace('Aaegypti_wt', 'Aedes_aegypti').replace('Aedessp', 'Aedes_sp')
    name = name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('//', '/').replace('__', '_').replace('.', '').replace(',', '')
    name = re.sub('^[\/\_\-]', '', name)
    try: # ID must start with letter
        name = 'V' + str(int(name))  # Deal with numbers (maybe python3)? V6845 ID requirement
    except:
        pass
    name = replace_strain_name(name, fix_whole_name)    # Before and after local processing? (try with and without, see if different)
    return name, original_name

# Load data
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

print(zika_fasta[0])

print("\n\nstrain: ",zika_fasta[0][0]['strain'])

print("fix_name output:", fix_name(zika_fasta[0][0]['strain']))