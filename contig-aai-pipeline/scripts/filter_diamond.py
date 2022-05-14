with open(snakemake.input[0]) as fin, open(snakemake.output[0], "w") as fout:
    current_qseqid = None
    for line in fin:
        qseqid, sseqid, _, qcovhsp, scovhsp, _ = line.split()
        qcovhsp, scovhsp = float(qcovhsp), float(scovhsp)
        qseqgen, sseqgen = qseqid.rsplit("_", 1)[0], sseqid.rsplit("_", 1)[0]
        if (
            (qseqgen != sseqgen)
            and (qcovhsp >= snakemake.params["min_cov"])
            and (scovhsp >= snakemake.params["min_cov"])
        ):
            if qseqid != current_qseqid:
                current_qseqid = qseqid
                current_qseqid_matches = set()
            if sseqgen not in current_qseqid_matches:
                current_qseqid_matches.add(sseqgen)
                fout.write(line)
