import Bio
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Align.Applications import ClustalwCommandline
import pandas as pd
from Bio import SeqIO
import config as cg
from collections import defaultdict
import os
import os.path
from Bio import Entrez

Entrez.email = "forexample@ucdconnect.ie"

# This two csv files were extracted from the supplementary material of this paper(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7011714/#R44)
# They are the reference sequences used for phylogenetic analysis (in cluding MSA) in RIVM
taxonomy = ['./new_data/rdrp_taxonomy.csv','./new_data/vp1_taxonomy.csv']

df = pd.read_csv(taxonomy[0], sep=",",skiprows = 1)

#Downloading all fastas for RdRp analysis
files = {}
for j in range(df.shape[0]):
    number = df.iloc[j]['GenBank accesion number']
    group = df.iloc[j]['New P-type'].split('.')[0]
    print(number)
    # Genbank doesn't have it, so skip it
    if number == 'F529737':
        continue
    # write fastas into respective groups
    if group not in files:
        files[group] = open('./new_data/reference_sequences_RdRp/'+group+'.fasta', 'w') 
    handle = Entrez.efetch(db="sequences", id = number, rettype="fasta", retmode="text")
    record = handle.read()
    files[group].write(record)
for file in files.values():
    file.close()

#Downloading all fastas for VP1 analysis
files = {}
df = pd.read_csv(taxonomy[1], sep=",",skiprows = 1)
for j in range(df.shape[0]):
    number = df.iloc[j]['GenBank accesion number']
    group = df.iloc[j]['New genotype'].split('.')[0]
    print(number)
    if number == 'F529737':
        continue
    if group not in files:
        files[group] = open('./new_data/reference_sequences_VP1/'+group+'.fasta', 'w') 
    handle = Entrez.efetch(db="sequences", id = number, rettype="fasta", retmode="text")
    record = handle.read()
    files[group].write(record)
for file in files.values():
    file.close()