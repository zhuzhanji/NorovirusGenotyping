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

# This is the accession id of the reference sequences for BLAST
# I collected them by analyzing the results of RIVM
numbers = {
    'GI':['NC_001959', 'NC_039897.1', 'NC_044853','NC_044854.1', 'NC_044856.1'],
    'GII':['NC_029646', 'NC_039475.1', 'NC_039476.1', 'NC_039477.1','NC_040876', 'NC_044045.1','NC_044046.1', 'NC_044932'], 
    'GIII':['NC_029645'],
    'GIV':['NC_029647', 'NC_044855.1'],
    'GV':['NC_008311', 'MW174170.1'],
    'GVI':['NC_044047','MW662289.1','MN908340.1','MW945229.1'],
    'GVII':['FJ692500', 'OL757872.1'],
    'GVIII':['AB985418.2'],
    'GIX':['OR050586.1','OR050585.1','OR050584.1'],
    'GX':['KJ790198.1','MF373609.1'],
    'Caliciviridae_Sapovirus_GI':['NC_006269'],
    'Caliciviridae_Sapovirus_GII':['NC_006554'],
    'Caliciviridae_Sapovirus_GIII':['NC_000940'],
    'Caliciviridae_Sapovirus_GIV':['DQ366346'],
    'Caliciviridae_Sapovirus_GV':['NC_027026'],
    'Caliciviridae_Sapovirus_GVI':['KJ508818'],
    'Caliciviridae_Sapovirus_GVII':['MF766258'],
    'Caliciviridae_Sapovirus_GVIII':['LC215895'],
    'Caliciviridae_Sapovirus_GX':['LC215896'],
    'Caliciviridae_Sapovirus_GXI':['LC215899'],
    'Caliciviridae_Sapovirus_GXII':['KX000385'],
    'Caliciviridae_Recovirus_unclassified':['EU391643'],
    'Caliciviridae_Chicken_calicivirus_unclassified':['HQ010042'],
    'Caliciviridae_Vesivirus_Feline_calicivirus':['NC_001481'],
    'Caliciviridae_Nacovirus_Turkey_calicivirus':['JQ347522'],
    'Caliciviridae_Lagovirus_Rabbit_hemorrhagic_disease_virus':['NC_001543'],
    'Caliciviridae_Lagovirus_unclassified':['NC_011704'],
    'Caliciviridae_unclassified_Bovine_calicivirus':['KT119483'],
    'Caliciviridae_Vesivirus_Vesicular_exanthema_of_swine_virus':['NC_002551'],
    'Caliciviridae_Lagovirus_European_brown_hare_syndrome_virus':['NC_002615'],
    'Caliciviridae_Vesivirus_Canine_calicivirus':['NC_004542'],
    'Caliciviridae_Nebovirus_unclassified':['NC_006875'],
    'Caliciviridae_Vesivirus_Rabbit_vesivirus':['NC_008580'],
    'Caliciviridae_Vesivirus_Steller_sea_lion_vesivirus':['NC_011050'],
    'Caliciviridae_Valovirus_unclassified':['NC_012699'],
    'Caliciviridae_Sapovirus_Bat_sapovirus':['NC_017936'],
    'Caliciviridae_Salmonid_calicivirus_unclassified':['NC_024031'],
    'Caliciviridae_Nacovirus_Goose_calicivirus':['NC_024078'],
    'Caliciviridae_Vesivirus_Ferret_badger_vesivirus':['NC_027122'],
    'Caliciviridae_Vesivirus_San_Miguel_sea_lion_virus':['U15301']
}

wholef = open('./new_data/database/refseq.fasta', 'w') 
for group, nums in numbers.items():
    for id in nums:
        handle = Entrez.efetch(db="sequences", id = id, rettype="fasta", retmode="text")
        record = handle.read()
        record.replace('>', '')
        wholef.write('>'+group+'|')
        wholef.write(record)
        #print(id)
wholef.close()
