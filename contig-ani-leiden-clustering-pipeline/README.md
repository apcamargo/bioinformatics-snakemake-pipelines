# contig-ani-leiden-clustering-pipeline

Computes *contig-wise* ANI from an input FASTA file, clusters the sequences, and writes cluster representatives.

```
# Run the ANI-based dereplication with default parameters
snakemake --config input=genomes.fasta --cores 8 -s contig-ani-leiden-clustering-pipeline.smk

# Change 'output', 'avg_ani', 'min_cov', 'blast_max_evalue', and 'blast_threads'
snakemake --config input=genomes.fasta \
                   output=representatives.fasta \
                   avg_ani=0.9 \
                   min_cov=0.5 \
                   blast_max_evalue=1e-10 \
                   blast_threads=16 \
                   --cores 16 -s contig-ani-leiden-clustering-pipeline.smk
```

## Dependencies
- Snakemake
- numpy
- igraph
- BLAST
- seqkit
