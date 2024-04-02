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


# This is the default prefix for all the fasta files
PATHS = {
    'GI':["./new_data/reference_sequences_RdRp/GI", "./new_data/reference_sequences_VP1/GI"],
    'GII':["./new_data/reference_sequences_RdRp/GII", "./new_data/reference_sequences_VP1/GII"],
    'GIII':["./new_data/reference_sequences_RdRp/GIII", "./new_data/reference_sequences_VP1/GIII"],
    'GIV':["./new_data/reference_sequences_RdRp/GIV", "./new_data/reference_sequences_VP1/GIV"],
    'GV':["./new_data/reference_sequences_RdRp/GV", "./new_data/reference_sequences_VP1/GV"],
    'GVI':["./new_data/reference_sequences_RdRp/GVI", "./new_data/reference_sequences_VP1/GVI"],
    'GVII':["./new_data/reference_sequences_RdRp/GVII", "./new_data/reference_sequences_VP1/GVII"],
    'GVIII':["", "./new_data/reference_sequences_VP1/GVIII"],
    'GNA1':["./new_data/reference_sequences_RdRp/GNA1", "./new_data/reference_sequences_VP1/GNA1"],
    'GNA2':["./new_data/reference_sequences_RdRp/GNA2", "./new_data/reference_sequences_VP1/GNA2"],
    'GIX':["", "./new_data/reference_sequences_VP1/GIX"],
    'GX':["./new_data/reference_sequences_RdRp/GX", "./new_data/reference_sequences_VP1/GX"],
}

# read fasta sequences into dictionary {id:fasta}
def ReadFastaFileAsDict(inpath: str) -> defaultdict(str):
    ori_seq = defaultdict(str)
    for seq_record in SeqIO.parse(inpath, "fasta"):
        id, sequence = seq_record.id, str(seq_record.seq)
        id = id.split('|')[-1]
        ori_seq[id] = sequence
    return ori_seq 

# For fastas under /rdrp directory, ORF1 coordinates will be extracted
# For fastas under /vp1 directory, ORF2 coordinates will be extracted
# Annotation files are written into '_annotation.csv' under the same directory of fasta files
for group, files in PATHS.items():
    for j in range(len(files)):
        orf1 = (j == 0)
        file = files[j] + '.fasta'
        if not os.path.isfile(file):
            continue
        ori_seq = ReadFastaFileAsDict(file)
        annotation = {}
        for id, seq in ori_seq.items():
            annotation[id] = getORF(id, orf1, not orf1)
            print(id, annotation[id])
        cols = ["rdrp_start", "rdrp_end"]
        if not orf1:
            cols = ["vp1_start", "vp1_end"]
        df = pd.DataFrame.from_dict(annotation, orient='index', columns=cols)
        df = df.reset_index().rename({'index':'id'},axis = 'columns')
        df.to_csv(files[j] + '_annotation.csv', sep=',', index=False, encoding='utf-8')  


# For fastas under /rdrp directory, ORF1 sequences will be extracted and written into .orf files
# For fastas under /vp1 directory, ORF2 sequences will be extracted and written into .orf files
for group, paths in PATHS.items():
    if os.path.isfile(paths[0] + '_annotation.csv'):
        # read the annotation file produced above
        df = pd.read_csv(paths[0] + '_annotation.csv', sep=",")
        file = open(paths[0] + '.orf', 'w')
        for seq_record in SeqIO.parse(paths[0]+'.fasta', "fasta"):
            id, sequence = seq_record.id.__str__(), seq_record.seq.__str__()
            # Write the orf1 fastas
            t1 = df[df['id'] == id].iloc[0]
            file.write('>'+id+'\n')
            file.write(sequence[max(t1['rdrp_start'] - 1, t1['rdrp_end'] - 762) :t1['rdrp_end']]+'\n')
        file.close()
        
    if os.path.isfile(paths[1] + '_annotation.csv'):
        df = pd.read_csv(paths[1] + '_annotation.csv', sep=",")
        # write orf2 fastas
        file = open(paths[1] + '.orf', 'w') 
        for seq_record in SeqIO.parse(paths[1]+'.fasta', "fasta"):
            id, sequence = seq_record.id.__str__(), seq_record.seq.__str__()            
            t1 = df[df['id'] == id].iloc[0]
            file.write('>'+id+'\n')
            file.write(sequence[t1['vp1_start'] - 1:t1['vp1_end']]+'\n')
        file.close()
