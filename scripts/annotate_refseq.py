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

# This functions retrives and parses the xml annotation file of a sequence (gbid) 
# gbid: genbank accension number, 
# orf1: if orf1 coordinates should be extracted, 
# orf2: if orf2 coordinates should be extracted,
# return value (list)
# {'orf1':[x,y], 'orf2':[x,y]}

def getORF(gbid:str, orf1: bool, orf2: bool):
    handle = Entrez.efetch(db="nucleotide", id = gbid, rettype="xml", retmode="xml")
    record = handle.read()
    tree = ET.fromstring(record)
    result = {}
    for child in tree.iter('GBFeature'):
        featurekey = child.find('GBFeature_key')
        if featurekey != None and featurekey.text == 'CDS':
            location = (child.find('GBFeature_location').text)
            for quals in child.iter('GBQualifier'):
                name = quals.find('GBQualifier_name').text
                if name == 'product':
                    val = quals.find('GBQualifier_value').text.lower()
                    if orf2:
                        if val.find('capsid') != -1 or val.find('vp1') != -1 :
                            result['orf2'] = location.split('..')
                    if orf1:
                        if val.find('polyprotein') != -1 or val.find('polymerase') != -1 or val.find('rdrp') != -1 or val.find('nonstructural protein') != -1 or val.find('orf1') != -1:
                            result['orf1'] = location.split('..')

    return result


# All the accession number for BLAST analysis
noronumbers = {
    'GI':['NC_001959', 'NC_039897.1', 'NC_044853','NC_044854.1', 'NC_044856.1'],
    'GII':['NC_029646', 'NC_039475.1', 'NC_039476.1', 'NC_039477.1','NC_040876', 'NC_044045.1','NC_044046.1', 'NC_044932'], 
    'GIII':['NC_029645'],
    'GIV':['NC_029647', 'NC_044855.1'],
    'GV':['NC_008311', 'MW174170.1'],
    'GVI':['NC_044047','MW662289.1','MN908340.1','MW945229.1'],
    'GVII':['FJ692500', 'OL757872.1'],
    'GVIII':['AB985418.2'],
    'GIX':['OR050586.1','OR050585.1','OR050584.1'],
    'GX':['KJ790198.1','MF373609.1']
}
annotation = {}
for group, nums in noronumbers.items():
    for id in nums:
        annotation[id] = [-1] * 4
        res = getORF(id, True, True)
        if 'orf1' in res:
            annotation[id][:2] = res['orf1']
        if 'orf2' in res:
            annotation[id][2:] = res['orf2']

# This annotation file will contain 5 columns: genbank id and coordinates
cols = ["rdrp_start", "rdrp_end", "vp1_start", "vp1_end"]
df = pd.DataFrame.from_dict(annotation, orient='index', columns=cols)
df = df.reset_index().rename({'index':'id'},axis = 'columns')
# Default location of this annotation file
df.to_csv("./new_data/database/refseq_annotation.csv", sep=',', index=False, encoding='utf-8') 