from pathlib import Path

ID = Path(config["input"]).stem

rule all:
    input: config.get("output", f"{ID}_aai.tsv")

rule prodigal:
    input: config["input"]
    output: f"{ID}_proteins.faa"
    threads: config.get("prodigal_threads", 8)
    script: "scripts/parallel_prodigal.py"

rule diamond:
    input: f"{ID}_proteins.faa"
    output: f"{ID}.dmd"
    params:
        max_eval = config.get("max_eval", 1e-3),
        min_seq_id = 100 * config.get("min_seq_id", 0.0)
    threads: config.get("diamond_threads", 8)
    shell:
        """
        diamond blastp --threads {threads} -e {params.max_eval} --sensitive \
        --id {params.min_seq_id} -q {input} -d {input} -o {output} \
        --max-target-seqs 25000 --outfmt 6 qseqid sseqid pident qcovhsp scovhsp length
        """

rule filter_diamond:
    input: f"{ID}.dmd"
    output: f"{ID}_filtered.dmd"
    params:
        min_cov = 100 * config.get("min_cov", 0)
    script: "scripts/filter_diamond.py"

rule get_rbh:
    input: f"{ID}_filtered.dmd"
    output: f"{ID}_rbh.dmd"
    script: "scripts/get_rbh.py"

rule aaicalc:
    input:
        FASTA=f"{ID}_proteins.faa",
        RBH=f"{ID}_rbh.dmd"
    output: config.get("output", f"{ID}_aai.tsv")
    params:
        min_frac_shared_genes = config.get("min_frac_shared_genes", 0.3),
        min_n_shared_genes = config.get("min_n_shared_genes", 0)
    script: "scripts/aaicalc.py"
