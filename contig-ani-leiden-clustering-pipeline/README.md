# contig-ani-leiden-clustering-pipeline

Computes *contig-wise* ANI from an input FASTA file, clusters the sequences, and writes cluster representatives.

```
# Run the ANI-based dereplication with default parameters
snakemake --config input=genomes.fasta --cores 8 -s contig-ani-leiden-clustering-pipeline.smk

# Change 'output', 'min_ani', 'min_cov', 'blast_max_evalue', and 'blast_threads'
snakemake --config input=genomes.fasta \
                   output=representatives.fasta \
                   leiden_resolution=0.7 \
                   min_ani=0.9 \
                   min_cov=0.5 \
                   blast_max_evalue=1e-10 \
                   blast_threads=16 \
                   --cores 16 -s contig-ani-leiden-clustering-pipeline.smk
```

The `leiden_resolution` parameter controls the granularity of the clusters. Higher values will result in smaller clusters that are highly interconnected.

This pipeline allows you to either set a specific value for the resolution parameter (default: `leiden_resolution=1.0`), or to let the algorithm choose a value for you (`leiden_resolution=auto`) by iterating over a range of values to minimize the difference between the average intra-cluster ANI and a target ANI (`avg_ani`). For the latter option, we recommend using `min_ani=0` to allow connections between sequences that have ANI lower than `avg_ani`.

```
snakemake --config input=genomes.fasta \
                   output=representatives.fasta \
                   leiden_resolution=auto \
                   avg_ani=0.95 \
                   min_ani=0.0 \
                   blast_threads=16 \
                   --cores 16 -s contig-ani-leiden-clustering-pipeline.smk
```

## Dependencies
- Snakemake
- numpy
- igraph
- BLAST
- seqkit
