## Part 1. Discoveries on the implementation of RIVM

### The reference dataset for BLAST analysis
In [A. Kroneman 2011], the authors wrote 'In the first step the query sequence is analyzed using BLAST against a reference set of whole genome sequences representing different taxonomic groups', but they didn't offer more details. in order to find out the reference dataset, a testset of over 5000 fastas of Caliciviridae family were uploaded to RIVM as queries. Result files were downloaded via the blue link at the bottom.

<img width="615" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/31eb7e55-b28b-4000-81ab-13e2de89cac5">

Here are the snapshots inside the result file. The column refseq contains all the reference sequences they use for BLAST analysis. It is a small reference dataset, one whole genome sequence (refseq) for one taxonomic group.

<img width="732" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/0efe39f0-86d8-470f-baf9-a127048c45dc">

<img width="726" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/70f30505-6298-454e-bb4e-cbbfe741951c">

### The reference dataset for MSA

Here are some snapshots of the report page. Here are what we can deduce from them:

1. They use PAUP* for phylogenetic tree construction
2. They use bootstrap

<img width="792" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/5accc785-a80e-4dfd-8e5a-08f2c5f1fb8a">

<img width="830" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/423382f5-4f29-487b-a84b-faf9655a5f3d">

<img width="379" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/8d1bdf04-b0db-47bc-8172-efc8fba882fa">


3. Here is the snapshot of the alignment. They align nucleotide sequences.

<img width="552" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/658a0af0-bd8f-45b2-972e-9044ba4a7437">


## Part2. Dependent packages
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


## Part3. Documentation of scripts

## data preparation

Original data and data produced are under 'new_data/', please refer to 'new_data/readme' for furthur knowledge.

### 1.Preparing Reference Sequences for BLAST analysis

- scripts/blast_reference_data.py

**function**

This standalone python script downloads all the reference sequences for BLAST analysis from GenBank. It then parses the GenBank XML files and extracts the coordinates of ORF1 and ORF2 for all the blast reference sequences. The results are written to "./new_data/database/refseq_annotation.csv". 

<img width="373" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/8c4a9a99-4404-4ce2-b4d1-38045073c9ae">


**Result**

All fastas are written into 'new_data/database/refseq_2.fasta'. The target location is hardcoded in 'DataPreparation.py'.

 <img width="723" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/2bfc3e6e-d3a7-4c6a-adb3-5f6cdaf9beb7">
 
The paser is **not 100% reliable**, some coordinates might contain leading or trailing '<','>', so manual correction is needed.




### 2.Preparing Reference Sequences for MSA from GenBank

- scripts/msa_reference_data.py

**Function**

This standalone python script downloads all the reference sequences for phylogenetic analysis from GenBank. 

**Result**

Fastas of the same genogroup will be written into 'group'.fasta, e.g. GII.fasta. 

<img width="108" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/6a804f09-9335-4c34-87f0-e18395c2f6dc">

Current locations are './new_data/reference_sequences_VP1/' for ORF2 reference sequences and './new_data/reference_sequences_RdRp/' for ORF1 reference sequences.


- scripts/extract_orfs.py

**Function**

This standalone python script retrieves and parses the GenBank XML files, extracts the coordinates of ORF1 and ORF2 for all the reference sequences for MSA. 

**Result**
The coordinates are written into '#GROUP_annotation.csv' ( e.g. GII_annotation.csv). 

<img width="156" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/c05a9068-0c1f-41fe-aca4-5a01717aef24">

ORF1 and ORF2 sequences will be extracted separately and saved as '#group.orf' under the same directory of fasta files.
All orf1 or rdrp region are written into rdrp/all.orf. All orf2 or rdrp region are written into vp1/all.orf. 

<img width="97" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/c7b6943e-48d5-4990-9a88-c98340ee6cdb">
<img width="686" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/ae393b49-1e31-485d-aa8c-e8cdffdcf5e7">


### 3.Translation of VP1 gene

**Function**

In DataPreparation.py, function TranslateORF2() translates all the ORF2 sequences into amino acids. 

**Result**

'All.pro' contains all the amino acids of ORF2 sequences of different genogroups. 

AA of different genogroups will be written into respective .pro files.

### 4.Making BLAST database
-scripts/makedb.sh
-scripts/makedb.py

**Function**

makedb.sh calls makedb.py. 

makedb.py takes two parameters: source fasta file, nucleotide blast database path.  

The default source fasta file is the refseq_2.fasta produced in step 1. 

The default name of the database is 'refseqdb_2', which is also under new_data/database.

**Result**

It makes a nt database from the source fasta file.

### 5.Pre-alignment of RdRp partial sequences and VP1 partial sequences

-clustalw.sh
-clustalw.py

**Function**

clustalw.sh calls clustalw.py. 

clustalw.py uses clustalw to align all GI orf1s into rdrp/GI.aln, all GI vp1/orf2s into rdrp/GII.aln, all GII orf1s into vp1/GII.aln, all GII orf2s into rdrp/GII.aln, orf1s of all genogroups into rdrp/all.orf, orf2s of all genogroups into vp1/all.orf.

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
