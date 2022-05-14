from pathlib import Path

ID = Path(config["input"]).stem

rule all:
    input: config.get("output", f"{ID}_dereplicated_sequences.fasta")

rule makeblastdb:
    input: config["input"]
    output: temp(directory(f"{ID}_blastdb"))
    shell: "makeblastdb -dbtype nucl -in {input} -out {output}/blastdb"

rule blastn:
    input:
        FASTA=config["input"],
        BLASTDB=f"{ID}_blastdb"
    output: f"{ID}_blast.tsv"
    params:
        min_blast_ident = 100 * config.get("min_blast_ident", 0.0)
    threads: config.get("blast_threads", 8)
    shell:
        """
        blastn -task megablast -max_target_seqs 25000 -perc_identity {params.min_blast_ident} \
        -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue qlen slen" \
        -num_threads {threads} -query {input.FASTA} -db {input.BLASTDB}/blastdb -out {output}
        """

rule anicalc:
    input: f"{ID}_blast.tsv"
    output: f"{ID}_ani.tsv"
    params:
        blast_max_evalue = config.get("blast_max_evalue", 1e-3)
    script: "scripts/anicalc.py"

rule aniclust:
    input:
        FASTA=config["input"],
        ANI=f"{ID}_ani.tsv"
    output: f"{ID}_clusters.tsv"
    params:
        avg_ani = config.get("avg_ani", 0.95),
        min_cov = config.get("min_cov", 0.85),
        seed = config.get("seed", 1953)
    script: "scripts/aniclust.py"

rule get_cluster_representatives:
    input:
        FASTA=config["input"],
        CLUSTERS=f"{ID}_clusters.tsv"
    output:
        REPRESENTATIVES_LIST=temp(f"{ID}_cluster_rep.txt"),
        REPRESENTATIVES_FASTA=config.get("output", f"{ID}_dereplicated_sequences.fasta")
    shell:
        """
        awk '{{print $1}}' {input.CLUSTERS} > {output.REPRESENTATIVES_LIST} && \
        seqkit grep -f {output.REPRESENTATIVES_LIST} {input.FASTA} > {output.REPRESENTATIVES_FASTA}
        """
