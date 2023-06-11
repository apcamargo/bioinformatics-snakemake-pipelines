import math
from collections import defaultdict, deque
from bidict import bidict


def szudzik_pair(*numbers):
    if len(numbers) < 2:
        raise ValueError("Szudzik pairing function needs at least 2 numbers as input")
    elif any((n < 0) or (not isinstance(n, int)) for n in numbers):
        raise ValueError("Szudzik pairing function maps only non-negative integers")
    numbers = deque(numbers)
    n1 = numbers.popleft()
    n2 = numbers.popleft()
    if n1 != max(n1, n2):
        mapping = math.pow(n2, 2) + n1
    else:
        mapping = math.pow(n1, 2) + n1 + n2
    mapping = int(mapping)
    if not numbers:
        return mapping
    else:
        numbers.appendleft(mapping)
        return szudzik_pair(*numbers)


def szudzik_unpair(number, n=2):
    if (number < 0) or (not isinstance(number, int)):
        raise ValueError("Szudzik unpairing function requires a non-negative integer")
    if number - math.pow(math.floor(math.sqrt(number)), 2) < math.floor(
        math.sqrt(number)
    ):
        n1 = number - math.pow(math.floor(math.sqrt(number)), 2)
        n2 = math.floor(math.sqrt(number))
    else:
        n1 = math.floor(math.sqrt(number))
        n2 = (
            number
            - math.pow(math.floor(math.sqrt(number)), 2)
            - math.floor(math.sqrt(number))
        )
    n1, n2 = int(n1), int(n2)
    if n > 2:
        return szudzik_unpair(n1, n - 1) + (n2,)
    else:
        return n1, n2


gene_count = defaultdict(int)
genomes = []
genes = []
with open(snakemake.input["FASTA"]) as fin:
    for line in fin:
        if line.startswith(">"):
            genome_id = line.strip().split()[0][1:].rsplit("_", 1)[0]
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
        qseqgen, sseqgen = qseqid.rsplit("_", 1)[0], sseqid.rsplit("_", 1)[0]
        gene_pair = tuple(sorted((genes[qseqid], genes[sseqid])))
        gene_pair = szudzik_pair(*gene_pair)
        genome_pair = tuple(sorted((genomes[qseqgen], genomes[sseqgen])))
        genome_pair = szudzik_pair(*genome_pair)
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
        genome_1, genome_2 = szudzik_unpair(genome_pair)
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
