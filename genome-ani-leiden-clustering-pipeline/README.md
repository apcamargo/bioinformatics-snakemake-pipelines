# genome-ani-leiden-clustering-pipeline

Computes *genome-wise* ANI from input FASTA files, clusters the genomes, and writes cluster representatives.

```
# Run the ANI-based dereplication with default parameters
snakemake --config input=genomes_directory --cores 8 -s genome-ani-leiden-clustering-pipeline.smk

# Change 'output', 'avg_ani', 'min_cov', 'blast_max_evalue', and 'blast_threads'
snakemake --config input=genomes_directory \
                   output=representative_genomes_directory \
                   avg_ani=0.9 \
                   min_cov=0.5 \
                   blast_max_evalue=1e-10 \
                   blast_threads=16 \
                   --cores 16 -s genome-ani-leiden-clustering-pipeline.smk
```

## Dependencies
- Snakemake
- BLAST
- numpy
- igraph
- csvtk
