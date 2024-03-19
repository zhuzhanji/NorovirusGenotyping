# A Commandline Norovirus Typing Tool
## 1. Introduction
Norovirus is a highly contagious virus that is a common cause of acute gastroenteritis, often referred to as the stomach flu or winter vomiting bug. It belongs to the family Caliciviridae and is a non-enveloped, single-stranded RNA virus. Norovirus is a significant public health concern due to its high prevalence, ability to cause outbreaks, and impact on individuals and communities. With no specific antiviral treatment, it causes an estimated 200,000 deaths per year, including 50,000 child deaths[^1]. Norovirus exhibits notable genetic diversity, with 10 genogroups, GI–GX, and >48 genotypes identified[^2]. Analyzing high throughput data of norovirus is crucial for gaining insights into its mechanisms of infection and transmission. However, existing taxonomic classification tools for norovirus are server-based, presenting challenges in both data throughput and user-friendliness. There is a pressing need for an efficient, locally-operated, command-line-based taxonomic classification tool for norovirus.

#### 1.1 Taxonomic classification of caliciviridae
#### 1.2 Discovery and genus classification of norovirus
#### 1.3 Norovirus Epidemiology
#### 1.4 History of viral taxonomy in general(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7157452/)
#### 1.5 History of norovirus taxonomy
##### 1.5.1 Methods applied and how they work
##### 1.5.2 Challenges within current taxonomy framework for high throughput sequencing
##### 1.5.3 Detailed Overview of current norovirus taxonomy tools


# Methods
305 nucleotide sequences for VP1 region phylogenetic analyses, and 232 nucleotide sequences for RdRp region phylogenetic analyses, representing the genetic diversity of all norovirus genogroups and genotypes were downloaded from GenBank, as described previously[^3]. (last download date: March 16 2024). 305 Complete ORF2 sequences (ranging from 530 to 580 amino acids in length) were extracted and translated from the former using Biopython[^4]. 232 partial nucleotide sequences (762 nucleotides) of the RdRp region at the 3′-end of ORF1 were extracted from the latter.



![sequences](https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/7c695662-7de8-477f-9346-4a30dc7d1865  "sequences" )

In addtion, 30 norovirus whole genome sequences from RefSeq dataset[^5] and GenBank were also downloaded.

<img width="600" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/bcc44256-9c9f-4d1c-a64a-97f0f459f067">


## Workflow

<img width="661" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/877673d4-159c-4077-9f54-d5d21de959ba">

<img width="710" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/7145f545-b53e-4589-a709-7cb72dbc681d">


# Results
## Genogrouping

## Genotyping

# Discussion


[^1]:World Health Organization. (2022). Norovirus. Retrieved from https://www.who.int/teams/immunization-vaccines-and-biologicals/diseases/norovirus.
[^2]:Chhabra P, de Graaf M, Parra GI, Chan MC, Green K, Martella V, et al. Updated classification of norovirus genogroups and genotypes. [Erratum in: J Gen Virol. 2020;101:893]. J Gen Virol. 2019;100:1393–406. 10.1099/jgv.0.001318.
[^3]:Chhabra, P. et al. Updated classification of norovirus genogroups and genotypes. J. Gen. Virol. 100, 1393–1406 (2019).
[^4]:Cock, P. J. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25, 1422–1423 (2009).
[^5]:Brister JR, Ako-Adjei D, Bao Y, Blinkova O. NCBI viral genomes resource. Nucleic Acids Res. 2015 Jan;43(Database issue):D571-7 PubMed PubMedCentral



