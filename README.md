# NorovirusGenotyping
Norovirus is a highly contagious virus that is a common cause of acute gastroenteritis, often referred to as the stomach flu or winter vomiting bug. It belongs to the family Caliciviridae and is a non-enveloped, single-stranded RNA virus. Norovirus is a significant public health concern due to its high prevalence, ability to cause outbreaks, and impact on individuals and communities. With no specific antiviral treatment, it causes an estimated 200,000 deaths per year, including 50,000 child deaths. Norovirus exhibits notable genetic diversity, with 10 genogroups, GI–GX, and >48 genotypes identified. Analyzing high throughput data of norovirus is crucial for gaining insights into its mechanisms of infection and transmission. However, existing taxonomic classification tools for norovirus are server-based, presenting challenges in both data throughput and user-friendliness. There is a pressing need for an efficient, locally-operated, command-line-based taxonomic classification tool for norovirus.

<img width="698" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/8340031f-169a-4a79-b11f-baaa2e185241">

<img width="709" alt="image" src="https://github.com/zhuzhanji/NorovirusGenotyping/assets/37281560/90cc15fc-9b48-4876-9b2c-f23858c8b735">

Problem:
1.What is the BLAST score on rivm.nl? It has different scale from the BLAST bitscore.

2.What if the top 1 hit only contains ORF1 or ORF2 region, but the query sequence contains both? Is it possible? (If this happens, MSA against ORF1 or ORF1 might not be conducted.) 

solution: Top N hit? How to decide N? What threshold? 

3.Some reference sequences for RdRp contain VP1 region but are not included as the reference sequences for VP1, as a result, their VP1 region are not annotated.

Solution: reference sequences for RdRp and VP1 should all be viruligned against ORF1 and ORF2 xml

4.PAUP* analysis. How to automatically assign a genotype? Parsing the NEXUS output file? 
