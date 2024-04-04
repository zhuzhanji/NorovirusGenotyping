## dependent packages
- blast
- clustalw
- python3
- 
After blast and clustalw are installed in your computer, you should check scripts/config.py.
The absolute paths of blastn, makeblastdb, clustalw2 are written in scripts/config.py, you should change them first to match the locations of these executive files.

<img width="719" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/af895850-6e48-4ba9-9850-ede72536c833">

You can also add blastn, makeblastdb, clustalw2 to the environment paths, then change "exe" list into this:
    "exe":
    {
        "makeblastdb":"makeblastdb",
        "blastn":"blastn",
        "clustal":"clustalw2"
    },



## data preparation

Original data and data produced are under 'new_data/', please refer to 'new_data/readme' for furthur knowledge.

### 1.Preparing Reference Sequences for BLAST analysis

- scripts/blast_reference_data.py

**function**

This standalone python script downloads all the reference sequences for BLAST analysis from GenBank. It then parses the GenBank XML files and extracts the coordinates of ORF1 and ORF2 for all the blast reference sequences. The results are written to "./new_data/database/refseq_annotation.csv". 

**Result**

All fastas are written into 'new_data/database/refseq_2.fasta'. The target location is hardcoded in 'DataPreparation.py'.
It should be noted that the paser is not 100% reliable, some coordinates might contain leading or trailing '<' or '>', so manual correction is needed.
 

### 2.Preparing Reference Sequences for MSA from GenBank

- scripts/msa_reference_data.py

**Function**

This standalone python script downloads all the reference sequences for phylogenetic analysis from GenBank. 

**Result**

Fastas of the same genogroup will be written into 'group'.fasta, e.g. GII.fasta. 

The current destination locations are './new_data/reference_sequences_VP1/' for ORF2 reference sequences and './new_data/reference_sequences_RdRp/' for ORF1 reference sequences.

- scripts/extract_orfs.py

**Function**

This standalone python script parses the GenBank XML files and gets the coordinates of ORF1 and ORF2 for all the reference sequences for MSA. 

**Result**
The coordinates are written into '#GROUP_annotation.csv' ( e.g. GII_annotation.csv). 
ORF1 and ORF2 sequences will be extracted separately and saved as '#group.orf' under the same directory of fasta files.


### 3.Translation of VP1 gene

**Function**

In DataPreparation.py, function TranslateORF2() translates all the ORF2 sequences into amino acids. 

**Result**

'All.pro' contains all the amino acids of ORF2 sequences of different genogroups. 

AA of different genogroups will also be written into different .pro files.

### 4.Making BLAST database
-scripts/makedb.sh
-scripts/makedb.py

**Function**

makedb.sh calls makedb.py. makedb.py takes two parameters: source fasta file, output database path.  

The default source fasta file is refseq_2.fasta produced in step 1. 

The default name of the database is 'refseqdb_2'. It's under new_data/database.

**Result**

It makes a nt database from the source fasta file.

### 5.Pre-alignment of RdRp partial sequences and VP1 partial sequences

-clustalw.sh
-clustalw.py

clustalw.sh calls clustalw.py. 

**Function**

clustalw.py uses clustalw to align all GI orf1s into rdrp/GI.aln, all GI vp1/orf2s into rdrp/GII.aln, all GII orf1 into vp1/GII.aln, all GII orf2 into rdrp/GII.aln, orf1 of all genogroups into vp1/all.orf.

**Result**

These files can be found under new_data/reference_sequences_RdRp/ and new_data/reference_sequences_VP1/

-------------
## Workflow

### 1.BLAST query and genogroup assignment

-blastn.py

**Function**

It conducts BLAST analysis for a fasta file (1-N fastas)

It accepts two parameters: query fasta file path (1-N fastas), result path (.csv)

e.g. python blastn.py './../new_data/test_norovirus.fasta' './../new_data/test_result.csv'

-workflow.py (under development)

Currently it can assign genogroup after BLAST analysis. It then checks if the matching segment meets the length and region criteria. If all criteria are met, we should align the query sequence against the pre-aligned clustal file. Then phylogenectic trees should be constructed from this alignment. A genotype then can be assigned based on the supporting value.
