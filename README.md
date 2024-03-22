# A Commandline Norovirus Typing Tool
## 1. Introduction


#### 1.1 Taxonomic classification of caliciviridae
The family Caliciviridae includes viruses with a linear positive-sense RNA genome of 6.4–8.5 kb with the non-structural and structural proteins encoded by different ORFs. Virions are non-enveloped particles ranging from 27–40 nm in diameter with an icosahedral symmetry.[^1] Norovirus is one of the eleven genera in family Caliciviridae. In addition, unclassified caliciviruses have been detected in many animals.[^2]

#### 1.2 Discovery and genus classification of norovirus
Noroviruses are non-enveloped, single-stranded RNA viruses, belonging to the genus Norovirus in the family Caliciviridae. Most noroviruses contain three open reading frames (ORFs), except for murine noroviruses, which contain a fourth ORF. ORF1 codes for the polymerase, ORF2 codes for the capsid protein(VP1) and ORF3 codes for minor structural protein VP2.[^3] Like most viruses, recombination in noroviruses greatly increases their genetic variation and affect their phylogenetic clusterings.[^4]



<img width="483" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/9cecac0f-97e9-40a0-a14a-ee0c0fc5d84e" title="Norovirus Genome[3]">

Noroviruses were first isolated and identified in 1972 using electron microscopy, and were originally named the "Norwalk agent" after the place where an outbreak of acute gastroenteritis occurred 4 years ago.[^5] 



[^1]:Vinje J, et al.; ICTV Report Consortium. 2019. ICTV virus taxonomy profile: caliciviridae. J Gen Virol. 100(11):1469–1470.
[^2]:Wille M, Eden JS, Shi M et al. Virus-virus interactions and host ecology are associated with RNA virome structure in wild birds. Mol. Ecol. 2018, DOI: 10.1111/mec.14918.
[^3]:COMPARE Europe (2015). Reference genomes - Norovirus. Available at: https://www.compare-europe.eu/library/reference-genomes/norovirus (Accessed: March 16, 2024).
[^4]:Worobey M, Holmes EC. Evolutionary aspects of recombination in RNA viruses. J Gen Virol. 1999;80:2535–43.

[^5]:Knowable Magazine (2017). Norovirus: The perfect pathogen. Available at: https://knowablemagazine.org/content/article/health-disease/2017/norovirus-perfect-pathogen (Accessed: March 16, 2024).


#### 1.3 Norovirus Epidemiology
Noroviruses is a highly contagious virus that is a common cause of acute gastroenteritis, often referred as the stomach flu or winter vomiting bug. It is a significant public health concern due to its high prevalence, ability to cause outbreaks, and impact on individuals and communities. With no specific antiviral treatment, it causes an estimated 200,000 deaths per year, including 50,000 child deaths[^4].

[^4]:World Health Organization. (2022). Norovirus. Retrieved from https://www.who.int/teams/immunization-vaccines-and-biologicals/diseases/norovirus.


#### 1.4 History of viral taxonomy
Virus taxonomy addresses the classification of viruses into categories called taxa and the nomenclature for taxa.[^5] Several systematic classification of virus have been proposed since the discovery of virus. The current official virus taxonomy is maintained by the International Committee on Taxonomy of Viruses (ICTV), while Baltimore Classification, which was published in 1971, is also used today, as informal highest ranks of virus classification, in parallel with the official taxonomy.[^5]  

[^5]:Kuhn JH. Virus taxonomy. In: Bamford DH and Zuckerman M (eds). Encyclopedia of Virology, 4th edn. Oxford: Academic Press; 2021. pp. 28–37.

#### 1.5 History of norovirus taxonomy
It was recognized in the mid-1990s that norovirus strains should be organized into different genetic groups. 

At first, noroviruses were divided into genogroups and genotypes based on partial RdRp sequences.[^6] As more sequences became available, genotyping methods using cut-off values for percentage pair-wise similarity started to be used for norovirus taxonomy.[^7][^8] But the high time complexity of this method makes it difficult to handle large numbers of sequences. In addition, as prototype norovirus strains are often quite different from more recent strains due to accumulation of mutations, it is difficult to assign types using pairwise similarity cut-offs. 

In 2013 the Norovirus Classification Working Group (NCWG) proposed a standardized nomenclature and typing system using phylogenetic clustering of the complete VP1 amino acid sequences.[^9]In addition, dual typing (ORF1-RdRp=P type, ORF2=genotype) was proposed to include diversity of the partial RdRp sequences to address the increasing norovirus diversity.[^9] Using 2×standard deviation(sd) criteria to group sequences into separate clusters, the number of genogroups of noroviruses were expanded to 10 (GI-GX) and the number of genotypes were expanded to 49. Based on nucleotide diversity in the RdRp region, noroviruses were divided into more than 60 P-types.[^10] 

It should be noted that, though dual typing is widely used around the world now, most norovirus sequences on Genbank have not been updated according to this nomenclature system and many of them have been misclassified. 

[^6]:Vinjé J, Koopmans MP. Molecular detection and epidemiology of small round-structured viruses in outbreaks of gastroenteritis in the Netherlands. J Infect Dis. 1996;174:610–615. doi: 10.1093/infdis/174.3.610.
[^7]:Vinje J,Hamidjaja RA,Sobsey MD. Development and application of acapsid VP1(region D) based reverse transcription PCR assay for genotyping of genogroup I and II noroviruses. J Virol Methods 2004;116(March(2)):109–17.
[^8]:Zheng DP, Ando T, Fankhauser RL, Beard RS, Glass RI, Monroe SS. Norovirus classification and proposed strain nomenclature. Virology 2006;346(March(2)):312–23.
[^9]:Kroneman A, Vega E, Vennema H, Vinjé J, White PA, et al. Proposal for a unified norovirus Nomenclature and genotyping. Arch Virol. 2013;158:2059–2068. doi: 10.1007/s00705-013-1708-5. 
[^10]:Chhabra, P. et al. Updated classification of norovirus genogroups and genotypes. J. Gen. Virol. 100, 1393–1406 (2019).



#### 1.6 Challenges within current taxonomy framework for high throughput sequencing
During the last decade, High Throughput Sequencing(HTS), also refered to as next generation sequencing, has greatly advanced our knowledge of viruses. However, it is important to consider the challenges HTS faces in norovirus taxonomy. RNA viruses like noroviruses are difficult to sequence and characterize using HTS. Firstly, they are short in genome length, approximately 7.5k. Secondly, noroviruses exhibit notable genetic diversity. Thirdly, they lack universally conserved markers and genome plasticity, which contributes challenges for common PCR-based approaches.[^11]

[^11]:Fitzpatrick, A.H.; Rupnik, A.; O’Shea, H.; Crispie, F.; Keaveney, S.; Cotter, P. High Throughput Sequencing for the Detection and Characterization of RNA Viruses. Front. Microbiol. 2021, 12, 190.

  
#### 1.7 Current norovirus taxonomy tools

Analyzing high throughput data of norovirus is crucial for gaining insights into its mechanisms of infection and transmission. Two gold-standard norovirus typing tools are available, the NoroNet typing tool[^12] and the human calicivirus typing tool (HuCaT)[^13].

The Norovirus Typing Tool (https://www.rivm.nl/mpf/norovirus/typingtool) uses a BLAST algorithm against a set of reference sequences followed by phylogenetic analysis to assign norovirus genotypes and P-types.[^12] 
[^12]:Kroneman A, Vennema H, Deforche K, v d Avoort H, Penaranda S, Oberste MS, et al., An automated genotyping tool for enteroviruses and noroviruses, J. Clin. Virol 51 (2011) 121–125. 

Hucat([HuCaT](https://calicivirustypingtool.cdc.gov/bctyping.html)) uses a set of curated reference sequences that are compared to query sequences using a k-mer (DNA substring) based algorithm. Outputs include alignments and phylogenetic trees of the 12 top matching reference sequences for each query.[^13]

[^13]:Tatusov RL, Chhabra P, Diez-Valcarce M, Barclay L, Cannon JL, Vinjé J. 2021. Human calicivirus typing tool: a web-based tool for genotyping human norovirus and sapovirus sequences. J Clin Virol 134:104718. https://doi.org/10.1016/j.jcv.2020.104718.


However, both tools are web-based, creating a break in the processing of sequencing data. There is a pressing need for an efficient, locally-operated, command-line-based taxonomic classification tool for norovirus.



# Methods
305 nucleotide sequences for VP1 region phylogenetic analyses, and 232 nucleotide sequences for RdRp region phylogenetic analyses, representing the genetic diversity of all norovirus genogroups and genotypes were downloaded from GenBank, as described previously[^12]. (last download date: March 16 2024). 305 Complete ORF2 sequences (ranging from 530 to 580 amino acids in length) were extracted and translated from the former using Biopython[Cock, P. J. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25, 1422–1423 (2009).]. 232 partial nucleotide sequences (762 nucleotides) of the RdRp region at the 3′-end of ORF1 were extracted from the latter.


![sequences](https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/7c695662-7de8-477f-9346-4a30dc7d1865  "sequences" )

In addtion, 30 norovirus whole genome sequences from RefSeq dataset[Brister JR, Ako-Adjei D, Bao Y, Blinkova O. NCBI viral genomes resource. Nucleic Acids Res. 2015 Jan;43(Database issue):D571-7 PubMed PubMedCentral] and GenBank were also downloaded.

<img width="600" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/bcc44256-9c9f-4d1c-a64a-97f0f459f067">


## Workflow

<img width="661" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/877673d4-159c-4077-9f54-d5d21de959ba">

<img width="710" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/7145f545-b53e-4589-a709-7cb72dbc681d">


# Results
## Genogrouping

## Genotyping

# Discussion





