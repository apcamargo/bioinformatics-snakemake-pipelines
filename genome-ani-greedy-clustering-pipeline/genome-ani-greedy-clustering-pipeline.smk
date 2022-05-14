import shutil
from pathlib import Path

ID = Path(config["input"])

rule all:
    input: config.get("output", f"{ID}_dereplicated")

rule concatenate_genomes:
    input: config["input"]
    output: f"{ID}_concatenated.fna",
    script: "scripts/concatenate_genomes.py"

rule makeblastdb:
    input: f"{ID}_concatenated.fna"
    output: temp(directory(f"{ID}_blastdb"))
    shell: "makeblastdb -dbtype nucl -in {input} -out {output}/blastdb"

rule blastn:
    input:
        FASTA=f"{ID}_concatenated.fna",
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

rule process_blast:
    input: f"{ID}_blast.tsv"
    output: temp(f"{ID}_blast_processed.tsv")
    script: "scripts/process_blast.py"

rule sort_blast:
    input: f"{ID}_blast_processed.tsv"
    output: temp(f"{ID}_blast_sorted.tsv")
    shell:
        """
        csvtk sort -t -k 1 -k 2 -o {output} {input}
        """

rule anicalc:
    input:
        FASTA=f"{ID}_concatenated.fna",
        SORTED_BLAST=f"{ID}_blast_sorted.tsv"
    output: f"{ID}_ani.tsv"
    params:
        blast_max_evalue = config.get("blast_max_evalue", 1e-3)
    script: "scripts/anicalc.py"

rule aniclust:
    input:
        FASTA=f"{ID}_concatenated.fna",
        ANI=f"{ID}_ani.tsv"
    output: f"{ID}_clusters.tsv"
    params:
        min_ani = config.get("min_ani", 0.95),
        min_cov = config.get("min_cov", 0.85)
    script: "scripts/aniclust.py"

rule get_cluster_representatives:
    input:
        INPUT_DIR=config["input"],
        CLUSTERS=f"{ID}_clusters.tsv"
    output:
        REPRESENTATIVES_DIR=directory(config.get("output", f"{ID}_dereplicated"))
    run:
        rep_list = []
        with open(input.CLUSTERS) as fin:
            for line in fin:
                rep_list.append(line.split()[0])
        in_dir = Path(input.INPUT_DIR)
        out_dir = Path(output.REPRESENTATIVES_DIR)
        out_dir.mkdir()
        for rep in rep_list:
            source = list(in_dir.glob(rep + "*"))[0]
            target = out_dir.joinpath(source.name)
            shutil.copy(source, target)
