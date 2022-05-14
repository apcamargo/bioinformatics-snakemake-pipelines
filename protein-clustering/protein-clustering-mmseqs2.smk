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

rule sequencedb:
    input: f"{ID}_precluster_rep_seq.fasta"
    output: temp(directory(f"{ID}_seqdb"))
    shell:
        """
        mkdir {output}
        mmseqs createdb {input} {output}/{output}
        """

rule cluster:
    input: f"{ID}_seqdb"
    output:
        temp(directory(f"{ID}_clustdb")),
        temp(directory(f"{ID}_clust_tmp"))
    params:
        sensitivity=config.get("sensitivity", 7.5),
        max_eval = config.get("max_eval", 1e-3),
        min_aln_cov=config.get("min_aln_cov", 0.0),
        cluster_min_seq_id = config.get("cluster_min_seq_id", 0.0)
    threads: config.get("cluster_threads", 16)
    shell:
        """
        mkdir {output[0]}
        mmseqs cluster --threads {threads} -s {params.sensitivity} \
        -e {params.max_eval} -c {params.min_aln_cov} --cov-mode 0 \
        --cluster-mode 0 --max-seqs 5000 --min-seq-id {params.cluster_min_seq_id} \
        --cluster-reassign 1 {input}/{input} {output[0]}/{output[0]} {output[1]}
        """

rule createtsv:
    input:
        f"{ID}_seqdb",
        f"{ID}_clustdb"
    output: f"{ID}_clusters.tsv"
    shell:
        """
        mmseqs createtsv --full-header 1 {input[0]}/{input[0]} \
        {input[0]}/{input[0]} {input[1]}/{input[1]} {output} && \
        sed -i.bak -e 's/^\"//' -e 's/\"\t\"/\t/' -e 's/\"$//' {output} && \
        rm {output}.bak
        """

checkpoint family_membership:
    input: f"{ID}_clusters.tsv"
    output: directory("family_membership")
    params:
        min_family_size=config.get("min_family_size", 5),
    run:
        shell("mkdir {output}")
        family_number = 1
        for line in shell("csvtk cut -H -t -f 1 {input} | uniq -c | sort -n -r", iterable=True):
            family_size = int(line.strip().split()[0])
            if family_size >= params.min_family_size:
                family_name = f"{ID}_{family_number:06}"
                family_representative = line.split(maxsplit=1)[1]
                shell("csvtk grep -H -t -f 1 -p '\"{family_representative}\"' {input} | "
                      "csvtk cut -H -t -f 2 > {output}/{family_name}.txt")
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
