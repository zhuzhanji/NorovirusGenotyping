paths = {
    # annotation files for blast ref seqs
    "annotation":
    {
        "GI":"./data/database/GI_annotation.csv",
        "GII":"GII_annotation.csv"
    },
    # virulign results
    "virulign":
    {
        "GI":["./data/reference_sequences_RdRp/GI_annotation.fasta", "./data/reference_sequences_VP1/GI_annotation.fasta"],
        "GII":["GII_annotation_rdrp.fasta", "GII_annotation_vp1.fasta"]
    },
    "fasta":
    {
        "GI":["./data/reference_sequences_RdRp/GI.fasta", "./data/reference_sequences_VP1/GI.fasta"],
        "GII":["GII.fasta", "GII.fasta"],
        "Combined": "combined.fasta"
    },
    "database": "",
    "exe":
    {
        "makeblastdb":"",
        "blastn":"",
        "clustal":""
    },
    # clustal alignment
    "alignment":
    {
        "GI":[],
        "GII":[]
    },
    "taxonomy":
    {
        "GI":["./data/reference_sequences_RdRp/GI_taxonomy.tsv",
        "./data/reference_sequences_VP1/GI_taxonomy.tsv"
        ],
        "GII":["GII_taxonomy.tsv",]
    }
}
