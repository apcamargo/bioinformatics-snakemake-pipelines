# genome-aai-pipeline

Computes *genome-wise* AAI values from input FASTA files using reciprocal best hits between proteins predicted with Prodigal.

```
# Run the AAI computation with default parameters
snakemake --config input=genomes_directory --cores 8 -s genome-aai-pipeline.smk

# Change 'output', 'min_aai', 'min_cov', 'diamond_max_evalue', 'min_frac_shared_genes', 'min_n_shared_genes', 'prodigal_threads', 'diamond_threads'
snakemake --config input=genomes_directory \
                   output=aai.tsv \
                   min_aai=0.4 \
                   min_cov=0.5 \
                   diamond_max_evalue=1e-10 \
                   min_frac_shared_genes=0.5 \
                   min_n_shared_genes=5 \
                   prodigal_threads=16 \
                   diamond_threads=16 \
                   --cores 16 -s genome-aai-pipeline.smk
```

## Dependencies
- Snakemake
- bidict
- DIAMOND
- Prodigal
