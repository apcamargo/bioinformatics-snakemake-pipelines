from pathlib import Path
from Bio import SeqIO

ID = Path(config["input"]).stem

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.family_membership.get(**wildcards).output[0]
    return expand("family_alignments/{family}.fasta", family=glob_wildcards(
        Path(checkpoint_output).joinpath("{family}.txt")).family)

rule all:
    input: aggregate_input

rule clean_fasta:
    input: config["input"]
    output: temp(f"{ID}_cleaned.fasta")
    run:
        f = open(str(output), "w")
        record_id_list = set()
        for record in SeqIO.parse(str(input), "fasta"):
            record_id = record.id.replace(",", "_").replace("\"", "_").replace("'", "_")
            if record_id in record_id_list:
                record_id = record_id + "_2"
            seq = str(record.seq).upper().replace("*", "")
            f.write(f">{record_id}\n{seq}\n")
            record_id_list.add(record_id)
        f.close()

rule precluster:
    input: f"{ID}_cleaned.fasta"
    output:
        temp(f"{ID}_precluster_rep_seq.fasta"),
        temp(f"{ID}_precluster_all_seqs.fasta"),
        temp(f"{ID}_precluster_cluster.tsv"),
        temp(directory(f"{ID}_preclust_tmp"))
    params:
        precluster_min_seq_id=config.get("precluster_min_seq_id", 0.99),
    threads: config.get("precluster_threads", 16)
    shell:
        """
        mmseqs easy-linclust --threads {threads} --kmer-per-seq 100 -c 1.0 \
        --cluster-mode 2 --cov-mode 1 --min-seq-id {params.precluster_min_seq_id} \
        {input} {ID}_precluster {output[3]}
        """

rule diamond:
    input: f"{ID}_precluster_rep_seq.fasta"
    output: f"{ID}.dmd"
    params:
        diamond_max_eval = config.get("diamond_max_eval", 1e-3),
        diamond_min_seq_id = 100 * config.get("diamond_min_seq_id", 0.0)
    threads: config.get("diamond_threads", 16)
    shell:
        """
        diamond blastp --threads {threads} -e {params.diamond_max_eval} --very-sensitive \
        --id {params.diamond_min_seq_id} -q {input} -d {input} -o {output} \
        --outfmt 6 qseqid sseqid qlen slen nident mismatch evalue
        """

rule process_diamond:
    input: f"{ID}.dmd"
    output: f"{ID}.abc"
    params:
        min_aln_cov=config.get("min_aln_cov", 0.0)
    shell:
        """
        awk '{{if (($5 + $6 >= {params.min_aln_cov} * $3) && \
        ($5 + $6 >= {params.min_aln_cov} * $4)) {{print $1, $2, $7}}}}' \
        {input} > {output}
        """

rule mcxload:
    input: f"{ID}.abc"
    output:
        f"{ID}.mci",
        f"{ID}.tab"
    shell:
        """
        mcxload -abc {input} --stream-mirror --stream-neg-log10 \
        -stream-tf 'ceil(300)' -o {output[0]} -write-tab {output[1]}
        """

rule mcl:
    input:
        f"{ID}.mci",
        f"{ID}.tab"
    output: f"{ID}.mcl"
    params:
        mcl_inflation = config.get("mcl_inflation", 1.5)
    threads: config.get("mcl_threads", 16)
    shell:
        """
        mcl {input[0]} -use-tab {input[1]} -o {output} -te {threads} -I {params.mcl_inflation}
        """

checkpoint family_membership:
    input: f"{ID}.mcl"
    output: directory("family_membership")
    params:
        min_family_size=config.get("min_family_size", 5),
    run:
        shell("mkdir {output}")
        family_number = 1
        with open(str(input)) as fin:
            for line in fin:
                family_size = len(line.split("\t"))
                if family_size >= params.min_family_size:
                    family_name = f"{ID}_{family_number:06}"
                    output_file = Path(str(output)).joinpath(family_name + ".txt")
                    with open(output_file, "w") as fout:
                        fout.write("\n".join(line.split()))
                    family_number += 1
                else:
                    break

rule extract_family_sequences:
    input:
        f"{ID}_cleaned.fasta",
        "family_membership/{family}.txt"
    output: "family_sequences/{family}.fasta"
    shell: "seqkit grep {input[0]} -f {input[1]} -w 0 > {output}"

rule mafft:
    input:
        "family_sequences/{family}.fasta"
    output: "family_alignments/{family}.fasta"
    threads: config.get("mafft_threads", 4)
    shell:
        """
        mafft --quiet --anysymbol --thread {threads} --auto {input} | \
        seqkit seq -w 0 > {output}
        """
