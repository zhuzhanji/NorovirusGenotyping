import Bio
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
import pandas as pd
from Bio import SeqIO
import config as cg
from collections import defaultdict

def MakeDatabase(inpath: str, outpath: str, exepath: str = cg.paths['exe']['makeblastdb']):
    cline = NcbimakeblastdbCommandline(exepath, dbtype="nucl", input_file = inpath, out = outpath)
    cline()
    print('Database file has been writen to ' + outpath)

def QueryNcl(querypath: str, outpath: str, dbpath : str = cg.paths['database'], exepath: str = cg.paths['exe']['blastn']):
    cline = NcbiblastnCommandline(exepath, query = querypath, db = dbpath, evalue=0.00001, out= outpath, outfmt=6)
    cline()

def GetTopHit(resultpath: str, hit: int):   
    results = pd.read_csv(resultpath, sep="\t", header=None)
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
    if hit != 0:
        tophit = results.iloc[:hit, :]
    else:
        tophit = results

    return tophit
    
def PostQuery(resultpath: str):   
    results = pd.read_csv(resultpath, sep="\t", header=None)
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


    tophit = results.iloc[0]
    subject_group_id, subject_start, subject_end = tophit['subject'], tophit['subject_start'], tophit['subject_end']
    group_id = subject_group_id.split('##')
    if len(group_id) == 1:
        print('It is not norovirus, it is', group_id[0])
        return group_id[0], -1
    else:
        print('It is norovirus ', group_id[0])
    
    group, subject_seq = group_id
    # filter
    if subject_end - subject_start < 100:
        print('Matching length is shorter than 100')
        return

    type_to_annotation = cg.paths['annotation']
    # read annotation csv into datafrme
    annotation = pd.read_csv(type_to_annotation[group], sep=",", encoding='utf-8')
    print('HERE', annotation.columns[0], subject_seq)
    info = annotation[annotation[annotation.columns[0]] == subject_seq]
    if info.size == 0:
        print('There is not annotation for ', subject_seq)
        
    info = info.iloc[0]
    if info["rdrp_end"] != -1 and min(info["rdrp_end"], subject_end) - max(info["rdrp_start"], subject_start) >= 100:
        print('Valid ORF1')
    else:
        print('Sequence does not overlap sufficiently (>100 nucleotides) with ORF1')

    if info["vp1_start"] != -1 and min(info["vp1_end"], subject_end) - max(info["vp1_start"], subject_start) >= 100:
        print('Valid ORF2')
    else:
        print('Sequence does not overlap sufficiently (>100 nucleotides) with ORF2')


def ReadFastaFileAsDict(inpath: str) -> defaultdict(str):
    ori_seq = defaultdict(str)
    for seq_record in SeqIO.parse(inpath, "fasta"):
        id, sequence = seq_record.id, str(seq_record.seq)
        ori_seq[id] = sequence
    return ori_seq 
  
def ReadVirulignSeq(inpath, skip = 1) -> defaultdict(str):
    viru_seq = defaultdict(str)
    i = 0
    for seq_record in SeqIO.parse(inpath, "fasta"):
        id, sequence = seq_record.id, str(seq_record.seq)
        i += 1
        if i == skip:
            continue
        viru_seq[id] = (sequence.replace('-', '').lower())
    return viru_seq 



def ExportAnnotation(ori_seq, viru_seq_rdrp, viru_seq_vp1, outpath):
     
    annotation = {}
    for id, seq in ori_seq.items():
        if viru_seq_rdrp and id in viru_seq_rdrp:
            queryin_short = './tmp.fasta'
            testout_short = './tmp_result.csv'
            
            file = open(queryin_short, 'w')
            file.write('>'+id+'\n')
            file.write(viru_seq_rdrp[id]+'\n')
            file.close()
    
            QueryNcl(querypath = queryin_short, outpath = testout_short)
            res = GetTopHit(testout_short, 0)
            res = res[res['subject'] == 'GII##' + id]
            if res.size == 0:
                print('can not find rdrp', id)
                continue
            res = res.iloc[0]

            if res['subject'].find(id) == -1:
                print('Error! Blast subject is ', res['subject'], ' instead of rdrp ', id)
                annotation[id] = [-1, -1]
            else:    
                annotation[id] = [res['subject_start'], res['subject_end']]
            
            print(id, len(viru_seq_rdrp[id]), res['subject_end'] - res['subject_start'])
            
            #pos = seq.find(viru_seq_rdrp[id][:100])
            #if pos != -1:
            #    annotation[id] = [pos, pos + len(seq) - 1]
            #else:
            #    annotation[id] = [-1, -1]
            #    print('invalid rdrp ', id)
        else:
            annotation[id] = [-1, -1]

        if viru_seq_vp1 and id in viru_seq_vp1:
            queryin_short = './tmp.fasta'
            testout_short = './tmp_result.csv'
            
            file = open(queryin_short, 'w')
            file.write('>'+id+'\n')
            file.write(viru_seq_vp1[id]+'\n')
            file.close()
    
            QueryNcl(querypath = queryin_short, outpath = testout_short)
            res = GetTopHit(testout_short, 0)
            res = res[res['subject'] == 'GII##' + id]
            if res.size == 0:
                print('can not find vp1', id)
                continue
            res = res.iloc[0]

            if res['subject'].find(id) == -1:
                print('Error! Blast subject is ', res['subject'], ' instead of vp1 ', id)
                annotation[id].extend([-1, -1])
            else:    
                annotation[id].extend([res['subject_start'], res['subject_end']])
              
            print(id, len(viru_seq_vp1[id]), res['subject_end'] - res['subject_start'])
        else:
            annotation[id].extend([-1, -1])
        
    cols = ["rdrp_start", "rdrp_end", "vp1_start", "vp1_end"]   
    df = pd.DataFrame.from_dict(annotation, orient='index', columns=cols)
    df = df.reset_index().rename({'index':'id'},axis = 'columns')   
    df.to_csv(outpath, sep=',', index=False, encoding='utf-8')    

 