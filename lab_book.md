0.All data used are under 'new_data/', please read new_data/readme for furthur knowledge.

1.Preparing Reference Sequences for BLAST analysis

- scripts/blast_reference_data.py

  This standalone python script downloads all the reference sequences for BLAST analysis from GenBank. All fastas are written into './new_data/database/refseq.fasta'. The target location can be changed in the script.

2.Preparing Reference Sequences for MSA from GenBank

- scripts/msa_reference_data.py

  This standalone python script downloads all the reference sequences for phylogenetic analysis from GenBank. Fastas of the same genogroup will be written into 'group'.fasta, e.g. GII.fasta. The current destination locations are './new_data/reference_sequences_VP1/' for ORF2 reference sequences and './new_data/reference_sequences_RdRp/' for ORF1 reference sequences.
  
3.Translation of VP1 gene

4.Making BLAST database

5.Pre-alignment of RdRp partial sequences and VP1 partial sequences

6.BLAST query and genogroup assignment

