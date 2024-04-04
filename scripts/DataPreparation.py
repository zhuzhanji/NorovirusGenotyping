
import config as cg
import utils as utl
import os.path
import Bio
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
Entrez.email = "yang.liu2@ucdconnect.ie"

# download reference sequences from GenBank
def DownloadFastaFromGenbak():  
    wholef = open('./../new_data/database/refseq_2.fasta', 'w') 
    for group, nums in cg.numbers.items():
        for id in nums:
            handle = Entrez.efetch(db="sequences", id = id, rettype="fasta", retmode="text")
            record.replace('>', '')
            wholef.write('<'+group+'|')
            wholef.write(record)
            #print(id)
    wholef.close()

# Extract ORF1/2 coordinates from GenBank XML 
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

# Annotate sequences for C-typing and P-typing separately, and save as '#Group#_annotation.csv'
def MakeAnnotation(): 
    for group, files in cg.paths['prefix'].items():
        for j in range(len(files)):
            orf1 = (j == 0)
            file = files[j] + '.fasta'
            if not os.path.isfile(file):
                continue
            ori_seq = utl.ReadFastaFileAsDict(file)
            annotation = {}
            for id, seq in ori_seq.items():
                if orf1:
                    annotation[id] = getORF(id, orf1, False)['orf1']
                else:
                    annotation[id] = getORF(id, False, True)['orf2']
                print(id, annotation[id])
            cols = ["rdrp_start", "rdrp_end"]
            if not orf1:
                cols = ["vp1_start", "vp1_end"]
            df = pd.DataFrame.from_dict(annotation, orient='index', columns=cols)
            df = df.reset_index().rename({'index':'id'},axis = 'columns')
            df.to_csv(files[j] + '_annotation.csv', sep=',', index=False, encoding='utf-8') 

#extract orf1 and orf2 from sequences and save as '#group#_.orf' files
def ExtractOrfs():
    for group, paths in cg.paths['prefix'].items():
        if os.path.isfile(paths[0] + '_annotation.csv'):
            df = pd.read_csv(paths[0] + '_annotation.csv', sep=",")
            file = open(paths[0] + '.orf', 'w')
            for seq_record in SeqIO.parse(paths[0]+'.fasta', "fasta"):
                id, sequence = seq_record.id.__str__(), seq_record.seq.__str__()
                t1 = df[df['id'] == id].iloc[0]
                file.write('>'+id+'\n')
                file.write(sequence[max(t1['rdrp_start'] - 1, t1['rdrp_end'] - 762) :t1['rdrp_end']]+'\n')
            file.close()
        
        if os.path.isfile(paths[1] + '_annotation.csv'):
            df = pd.read_csv(paths[1] + '_annotation.csv', sep=",")
            file = open(paths[1] + '.orf', 'w') 
            for seq_record in SeqIO.parse(paths[1]+'.fasta', "fasta"):
                id, sequence = seq_record.id.__str__(), seq_record.seq.__str__()            
                t1 = df[df['id'] == id].iloc[0]
                file.write('>'+id+'\n')
                file.write(sequence[t1['vp1_start'] - 1:t1['vp1_end']]+'\n')
            file.close()


def TranslateORF2():
    taxfile2 = pd.read_csv(cg.paths['taxonomy']['vp1'], sep=",",skiprows = 1)
    orf2pall = open("./../new_data/reference_sequences_VP1/all.pro", 'w')
    for group, paths in cg.paths['prefix'].items():
        if os.path.isfile(paths[1] + '.orf'):
            orf2p = open(paths[1] + '.pro', 'w')
            for seq_record in SeqIO.parse(paths[1] + '.orf', "fasta"):
                id, sequence = seq_record.id.__str__(), (seq_record)
                if len(sequence) % 3 == 1:
                    sequence = sequence + Seq('nn')
                elif len(sequence) % 3 == 2:
                    sequence = sequence + Seq('n')
                sequence = (sequence.translate().seq).__str__()   
                t1 = taxfile2[taxfile2['GenBank accesion number'] == id.split('.')[0]]             
                orf2p.write('>' + group + '|' + t1.iloc[0]['New genotype'] + '|' + id+'\n')
                orf2p.write(sequence+'\n')

                orf2pall.write('>' + group + '|' + t1.iloc[0]['New genotype'] + '|' + id+'\n')
                orf2pall.write(sequence+'\n')

            orf2p.close()
    orf2pall.close()


#TranslateORF2()