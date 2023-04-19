from collections import namedtuple

Hsp = namedtuple(
    "Hsp",
    ["qname", "tname", "pid", "len", "qcoords", "tcoords", "evalue", "qlen", "tlen"],
)


def parse_blast(handle):
    for line in handle:
        line = line.split()
        yield Hsp(
            qname=line[0],
            tname=line[1],
            pid=float(line[2]) / 100,
            len=float(line[3]),
            qcoords=sorted([int(line[4]), int(line[5])]),
            tcoords=sorted([int(line[6]), int(line[7])]),
            evalue=float(line[8]),
            qlen=float(line[9]),
            tlen=float(line[10]),
        )


def yield_alignment_blocks(handle):
    # init block with 1st record
    key, alns = None, None
    for aln in parse_blast(handle):
        key = (aln.qname, aln.tname)
        alns = [aln]
        break
    # loop over remaining records
    for aln in parse_blast(handle):
        # extend block
        if (aln.qname, aln.tname) == key:
            alns.append(aln)
        # yield block and start new one
        else:
            yield alns
            key = (aln.qname, aln.tname)
            alns = [aln]
    yield alns


def prune_alns(alns, max_evalue=snakemake.params.blast_max_evalue):
    # remove alignments with > max_evalue
    # discard alignments after the entire query length has been covered
    keep = []
    cur_aln = 0
    qry_len = alns[0].qlen
    for aln in alns:
        qcoords = aln.qcoords
        aln_len = max(qcoords) - min(qcoords) + 1
        if aln.evalue > max_evalue:
            continue
        if cur_aln >= qry_len or aln_len + cur_aln >= 1.10 * qry_len:
            break
        keep.append(aln)
        cur_aln += aln_len
    return keep


def compute_ani(alns):
    return round(sum(a.len * (a.pid) for a in alns) / sum(a.len for a in alns), 4)


def compute_cov(alns):
    # merge qcoords
    coords = sorted([a.qcoords for a in alns])
    nr_coords = [coords[0]]
    for start, stop in coords[1:]:
        # overlapping, update start coord
        if start <= (nr_coords[-1][1] + 1):
            nr_coords[-1][1] = max(nr_coords[-1][1], stop)
        # non-overlapping, append to list
        else:
            nr_coords.append([start, stop])
    # compute qry_cov
    alen = sum(stop - start + 1 for start, stop in nr_coords)
    qcov = round(alen / alns[0].qlen, 4)
    # merge tcoords
    coords = sorted([a.tcoords for a in alns])
    nr_coords = [coords[0]]
    for start, stop in coords[1:]:
        # overlapping, update start coord
        if start <= (nr_coords[-1][1] + 1):
            nr_coords[-1][1] = max(nr_coords[-1][1], stop)
        # non-overlapping, append to list
        else:
            nr_coords.append([start, stop])
    # compute qry_cov
    alen = sum(stop - start + 1 for start, stop in nr_coords)
    tcov = round(alen / alns[0].tlen, 4)
    return qcov, tcov


out = open(snakemake.output[0], "w")
out.write("\t".join(["contig_1", "contig_2", "num_alns", "ani", "qcov", "tcov"]) + "\n")
input = open(snakemake.input[0])
for alns in yield_alignment_blocks(input):
    alns = prune_alns(alns)
    if len(alns) == 0:
        continue
    qname, tname = alns[0].qname, alns[0].tname
    ani = compute_ani(alns)
    qcov, tcov = compute_cov(alns)
    row = [qname, tname, len(alns), ani, qcov, tcov]
    out.write("\t".join(str(_) for _ in row) + "\n")
input.close()
out.close()
