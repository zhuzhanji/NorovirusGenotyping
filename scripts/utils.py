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

def MakeDatabase(inpath: str, outpath: str, exepath: str = cg.paths['exe']['makeblastdb']):
    cline = NcbimakeblastdbCommandline(exepath, dbtype="nucl", input_file = inpath, out = outpath)
    cline()
    print('Database file has been writen to ' + outpath)

def QueryNcl(querypath: str, outpath: str, dbpath : str = cg.paths['database'], exepath: str = cg.paths['exe']['blastn']):
    tmp_out = './out'
    cline = NcbiblastnCommandline(exepath, query = querypath, db = dbpath, evalue=0.00001, out= tmp_out, outfmt=6, \
                        max_target_seqs = 1, task = 'blastn', max_hsps = 1)
    cline()

    results = pd.read_csv(tmp_out, sep="\t", header=None)
    if results.shape[0] == 0:
        print('Blast Query Failed')
        return
    # Define column headers
    headers = ['query', 'subject',
           'pc_identity', 'aln_length', 'mismatches', 'gaps_opened',
           'query_start', 'query_end', 'subject_start', 'subject_end',
           'e_value', 'bitscore']
    # Assign headers
    results.columns = headers

    results.insert(2, "BLAST", '')
    for i in range(len(results)):
        subject, e_value = results.loc[i, 'subject'], results.loc[i,'e_value']
        if e_value >= 1e-5:
            results.loc[i, 'BLAST'] = 'unassigned'
        else:
            results.loc[i, 'BLAST'] = subject.split('|')[0]

    results.to_csv(outpath, sep='\t', index=False, encoding='utf-8') 
    os.remove(tmp_out)



def ReadFastaFileAsDict(inpath: str) -> defaultdict(str):
    ori_seq = defaultdict(str)
    for seq_record in SeqIO.parse(inpath, "fasta"):
        id, sequence = seq_record.id, str(seq_record.seq)
        id = id.split('|')[-1]
        ori_seq[id] = sequence
    return ori_seq 
  


 # Clustalw msa
def Clustalw(sequences: str, output: str, exepath: str = cg.paths['exe']['clustal']):
    cline = ClustalwCommandline(exepath, infile=sequences, outfile = output)
    print(cline)
    cline()

# Clustalw msa, sequence against alignment
def ClustalwAln(query: str, alignment: str, output: str, exepath: str = cg.paths['exe']['clustal']):
    # clustalw2 -profile1=file1.aln -profile2=query.fasta -sequences
    cline = ClustalwCommandline(exepath, profile1=alignment, profile2=query, sequences = True, outfile = output, output = 'PHYLIP')    
    cline()

def GetGenotype(id: str, taxfile: str):
    tax1 = pd.read_csv(taxfile, sep="\t|;", header=None, engine='python')
    t1 = tax1[tax1[0] == id]
    if len(t1) != 0:
        return t1[t1.shape[1] - 1].values[0]
    return ""



