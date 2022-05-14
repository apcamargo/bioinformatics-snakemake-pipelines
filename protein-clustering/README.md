# snakemake-pipelines

## `protein-clustering-diamond-mcl.smk`

Clusters and aligns proteins from an input FASTA file using DIAMOND for pairwise sequence comparison, MCL for clustering, and MAFFT for multiple sequence alignment.

```
# Run the clustering with the deafult parameters
snakemake --config input=proteins.faa -j 16 -s protein-clustering-diamond-mcl.smk

# Change 'precluster_min_seq_id', 'diamond_min_seq_id', 'min_aln_cov', and 'mcl_inflation'
snakemake --config input=proteins.faa \
                   precluster_min_seq_id=0.95 \
                   diamond_min_seq_id=0.8 \
                   min_aln_cov=0.8 \
                   mcl_inflation=3 \
                   -j 16 -s protein-clustering-diamond-mcl.smk
```

## `protein-clustering-mmseqs2.smk`

Clusters and aligns proteins from an input FASTA file using the set cover algorithm from MMSeqs2 for clustering and MAFFT for multiple sequence alignment.

```
# Run the clustering with the deafult parameters
snakemake --config input=proteins.faa -j 16 -s protein-clustering-mmseqs2.smk

# Change 'precluster_min_seq_id', 'min_aln_cov', and 'min_family_size'
snakemake --config input=proteins.faa \
                   precluster_min_seq_id=0.95 \
                   min_aln_cov=0.8 \
                   min_family_size=1 \
                   -j 16 -s protein-clustering-mmseqs2.smk
```
