import math
from collections import defaultdict
from bidict import bidict


def cantor_pair(k1, k2):
    z = int(0.5 * (k1 + k2) * (k1 + k2 + 1) + k2)
    return z


def cantor_depair(z):
    w = math.floor((math.sqrt(8 * z + 1) - 1) / 2)
    t = (w ** 2 + w) / 2
    y = int(z - t)
    x = int(w - y)
    return x, y


gene_count = defaultdict(int)
genomes = []
genes = []
with open(snakemake.input["FASTA"]) as fin:
    for line in fin:
        if line.startswith(">"):
            genome_id = line.strip().split()[0][1:].rsplit("_", 2)[0]
            gene_count[genome_id] += 1
            genomes.append(genome_id)
            genes.append(line.strip().split()[0][1:])
genomes = bidict({j: i for i, j in enumerate(sorted(set(genomes)))})
genes = bidict({j: i for i, j in enumerate(sorted(set(genes)))})

genome_pairs_genes = defaultdict(list)
gene_pair_info = defaultdict(lambda: [0, 0, 0])
with open(snakemake.input["RBH"]) as fin:
    for line in fin:
        qseqid, sseqid, pident, _, _, length = line.split()
        pident, length = float(pident), int(length)
        qseqgen, sseqgen = qseqid.rsplit("_", 2)[0], sseqid.rsplit("_", 2)[0]
        gene_pair = tuple(sorted((genes[qseqid], genes[sseqid])))
        gene_pair = cantor_pair(*gene_pair)
        genome_pair = tuple(sorted((genomes[qseqgen], genomes[sseqgen])))
        genome_pair = cantor_pair(*genome_pair)
        genome_pairs_genes[genome_pair].append(gene_pair)
        gene_pair_info[gene_pair][0] += pident
        gene_pair_info[gene_pair][1] += length
        gene_pair_info[gene_pair][2] += 1

with open(snakemake.output[0], "w") as fout:
    fout.write(
        "genome_1\t"
        "genome_2\t"
        "n_genes_genome_1\t"
        "n_genes_genome_2\t"
        "n_shared_genes\t"
        "aai\n"
    )
    for genome_pair, gene_pairs in genome_pairs_genes.items():
        genome_1, genome_2 = cantor_depair(genome_pair)
        genome_1 = genomes.inv[genome_1]
        genome_2 = genomes.inv[genome_2]
        min_n_genes = min(gene_count[genome_1], gene_count[genome_2])
        gene_pairs = set(gene_pairs)
        n_shared_genes = len(gene_pairs)
        if (
            n_shared_genes / min_n_genes >= snakemake.params["min_frac_shared_genes"]
        ) and (n_shared_genes >= snakemake.params["min_n_shared_genes"]):
            pair_aai = 0
            pair_total_length = 0
            for gene_pair in gene_pairs:
                pair_aai += (
                    gene_pair_info[gene_pair][0] / gene_pair_info[gene_pair][2]
                ) * (gene_pair_info[gene_pair][1] / gene_pair_info[gene_pair][2])
                pair_total_length += (
                    gene_pair_info[gene_pair][1] / gene_pair_info[gene_pair][2]
                )
            pair_aai /= pair_total_length
            pair_aai /= 100
            fout.write(
                f"{genome_1}\t"
                f"{genome_2}\t"
                f"{gene_count[genome_1]}\t"
                f"{gene_count[genome_2]}\t"
                f"{n_shared_genes}\t"
                f"{pair_aai:.4f}\n"
            )
