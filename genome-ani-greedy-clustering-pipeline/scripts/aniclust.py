from collections import defaultdict


def parse_seqs(path):
    with open(path) as handle:
        seq_id = next(handle).split()[0][1:]
        seq = ""
        for line in handle:
            if line[0] == ">":
                yield seq_id, seq
                seq_id = line.split()[0][1:]
                seq = ""
            else:
                seq += line.rstrip()
        yield seq_id, seq


# list seqs, sorted by length
seqs = defaultdict(int)
for seq_id, seq in parse_seqs(snakemake.input.FASTA):
    seq_id = seq_id.rsplit("_", 1)[0]
    # get the sequence length minus the number of Ns
    seqs[seq_id] += len(seq) - seq.upper().count("N")
seqs = [x[0] for x in sorted(seqs.items(), key=lambda x: x[1], reverse=True)]

# store edges
num_edges = 0
edges = dict([(x, []) for x in seqs])
with open(snakemake.input.ANI) as handle:
    for line in handle:
        qname, tname, num_alns, ani, qcov, tcov = line.split()
        if qname == tname:
            continue
        elif qname not in edges or tname not in edges:
            continue
        elif (
            float(tcov) < snakemake.params.min_cov
            or float(ani) < snakemake.params.min_ani
        ):
            continue
        edges[qname].append(tname)
        num_edges += 1

# cluster
clust_to_seqs = {}
seq_to_clust = {}
# loop over seqs in sorted order
for seq_id in seqs:
    # seq already been assigned; cant be centroid
    if seq_id in seq_to_clust:
        continue
    # add self to cluster
    clust_to_seqs[seq_id] = [seq_id]
    seq_to_clust[seq_id] = seq_id
    # update with cluster members
    for mem_id in edges[seq_id]:
        if mem_id not in seq_to_clust:
            clust_to_seqs[seq_id].append(mem_id)
            seq_to_clust[mem_id] = seq_id

# write output
out = open(snakemake.output[0], "w")
for seq_id, mem_ids in clust_to_seqs.items():
    out.write(seq_id + "\t" + ",".join(mem_ids) + "\n")
