from collections import defaultdict

with open(snakemake.input[0]) as fin:
    best_hits = defaultdict(dict)
    for line in fin:
        qseqid, sseqid, _, _, _, _ = line.split()
        qseqgen, sseqgen = qseqid.rsplit("_", 1)[0], sseqid.rsplit("_", 1)[0]
        best_hits[qseqid][sseqgen] = sseqid

with open(snakemake.input[0]) as fin, open(snakemake.output[0], "w") as fout:
    for line in fin:
        qseqid, sseqid, _, _, _, _ = line.split()
        qseqgen, sseqgen = qseqid.rsplit("_", 1)[0], sseqid.rsplit("_", 1)[0]
        if (qseqgen in best_hits[sseqid]) and (best_hits[sseqid][qseqgen] == qseqid):
            fout.write(line)
