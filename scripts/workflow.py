import utils as utl
import pandas as pd
import config as cg

def WorkFlow(querypath: str):
    # step1, Blast query against large database, assign genogroup
    utl.QueryNcl(querypath, './tmp.csv')
    results = pd.read_csv('./tmp.csv', sep="\t")
    if results.shape[0] == 0:
        print('Blast Query Failed')
        return

    tophit = results.iloc[0]
    subject_group_id, subject_start, subject_end = tophit['subject'], tophit['subject_start'], tophit['subject_end']
    group = tophit['BLAST']
    norogroup = ['GI','GII','GIII','GIV','GV','GVI','GVII','GVIII','GIX','GX']
    if group == 'unassigned' or group not in norogroup:
        return

    annotation = pd.read_csv(cg.paths['annotation'], sep=",")
    line = annotation[annotation['id'] == subject_group_id.split('|')[-1]].iloc[0] 

    if min(subject_end, line['rdrp_end']) - max(subject_start, line['rdrp_start']) < 100:
        print(subject_end, subject_start)
        print(line['rdrp_end'], line['rdrp_start'])
        print('Sequence does not overlap sufficiently (>100 nucleotides) with ORF1')
        
    else:
        #MSA aginast ORF1.aln
        # step 3.1
        #phylogenetic analysis
        #utl.ClustalwAln(querypath, "./../new_data/reference_sequences_RdRp/GI.aln", './EU085529_aln')
        print('ofr1 msa')

    if min(subject_end, line['vp1_end']) - max(subject_start, line['vp1_start']) < 100:
        print('Sequence does not overlap sufficiently (>100 nucleotides) with ORF2')
        print(subject_end, subject_start)
        print(line['rdrp_end'], line['rdrp_start'])
    else:
        #MSA aginast ORF2.aln
        # step 3.1
        #phylogenetic analysis
        print('orf2 msa')
        print(subject_end, subject_start)
        print(line['rdrp_end'], line['rdrp_start'])
        #utl.ClustalwAln(querypath, "./../new_data/reference_sequences_VP1/GI.aln", './EU085529_vp1_aln')
  

#EU085529.1
WorkFlow('./../testdata/AB089882.fasta')
