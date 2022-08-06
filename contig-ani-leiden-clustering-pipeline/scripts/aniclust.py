from collections import defaultdict

import itertools
import random
import igraph as ig
import numpy as np


class Sequence:
    def __init__(self, header, seq):
        self._header = header
        self._seq = seq

    @property
    def header(self):
        return self._header

    @property
    def id(self):
        return self._header.split()[0]

    @property
    def seq(self):
        return self._seq

    def count(self, substring):
        return self.seq.count(substring)

    def __len__(self):
        return len(self.seq)


def read_fasta(filepath):
    with open(filepath) as fin:
        last = None
        while True:
            if not last:
                for l in fin:
                    if l[0] == ">":
                        last = l[:-1]
                        break
            if not last:
                break
            name, seqs, last = last[1:], [], None
            for l in fin:
                if l[0] == ">":
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            seqs = "".join(seqs)
            if len(seqs):
                yield Sequence(name, seqs)
            if not last:
                break


def build_graph(
    ani_file,
    fasta_file,
    min_ani=snakemake.params.min_ani,
    min_cov=snakemake.params.min_cov,
):
    # Read the input FASTA file to get the names of all the sequences (nodes)
    nodes = []
    lengths = []
    for seq in read_fasta(fasta_file):
        nodes.append(seq.id)
        lengths.append(len(seq) - seq.seq.upper().count("N"))
    nodes_set = set(nodes)
    # Initiate the lists that will store the edges between two sequences and
    # their respective weights (ANI × coverage) and ANI
    edges = []
    weights = defaultdict(list)
    anis = defaultdict(list)
    with open(ani_file) as fin:
        # Skip header
        next(fin)
        for line in fin:
            (
                seq_1,
                seq_2,
                _,
                ani,
                qcov,
                tcov,
            ) = line.strip().split("\t")
            ani, qcov, tcov = (
                float(ani),
                float(qcov),
                float(tcov),
            )
            # If both the sequences in the ANI pair were in the input FASTA file,
            # their ANI is equal to or greater than `min_ani`, and either`qcov`
            # or `tcov` is equal to or greater than `min_cov`, store their connection
            if (
                ((seq_1 in nodes_set) and (seq_2 in nodes_set))
                and (ani >= min_ani)
                and ((qcov >= min_cov) or (tcov >= min_cov))
            ):
                pair = tuple(sorted([seq_1, seq_2]))
                edges.append(pair)
                anis[pair].append(ani)
                weight = ani * max([qcov, tcov])
                weights[pair].append(weight)
    # Take the mean ANI and weight for each pair of sequences
    anis = [np.mean(anis[pair]) for pair in edges]
    weights = [np.mean(weights[pair]) for pair in edges]
    # Create a graph from the node list and add edges weighted by the ANI between
    # the sequences being connected
    graph = ig.Graph()
    graph.add_vertices(list(nodes), attributes={"length": lengths})
    graph.add_edges(edges, attributes={"weight": weights, "ani": weights})
    return graph


def pick_resolution(
    graph,
    target_avg_ani=snakemake.params.avg_ani,
    steps=101,
    seed=snakemake.params.seed,
):
    """
    Given a graph (`graph`) and a target average within-cluster ANI (`target_avg_ani`),
    get the resolution parameter that will provide the closest within-community average
    edge weight using the Leiden clustering algorithm.
    """
    # The `last_res` variable will store the resolution (`res`) value of the previous
    # iteration, while the `last_avg_weight` will store the difference between the average
    # AAI of the previous iteration and the target average ANI (`target_avg_ani`)
    last_res = None
    last_avg_weight = None
    # Iterate through resolution parameters (`res`)
    for res in np.linspace(2, 0, steps):
        random.seed(seed)
        # Find the communities using the current resolution parameter
        communities = graph.community_leiden(weights="weight", resolution_parameter=res)
        # Initiate the list that will store the edges' ANI
        ani_list = []
        # Iterate through the communities subgraphs ordered by length (largest → smallest)
        for sg in reversed(np.argsort(communities.sizes())):
            members = [i.attributes()["name"] for i in communities.subgraph(sg).vs]
            # If the community has a single member, break the loop
            if len(members) < 2:
                break
            # Add all the pairwise ANI value to the `ani_list` list
            for i, j in itertools.combinations(members, 2):
                try:
                    edge = graph.get_eid(i, j)
                    edge = graph.es[edge]
                    ani_list.append(edge.attributes()["ani"])
                # In case two nodes are not connected, skip the pair
                except Exception:
                    continue
        # Compute the average value of all the edges in the iteration
        current_avg_ani = np.mean(ani_list) if len(ani_list) else 1
        # If this is the first iteration (that is, `res == 1`)
        if res == 2:
            last_res = res
            last_avg_weight = current_avg_ani
        # If the difference between the curent iteration average ANI and the target ANI
        # is less than the difference of the previous iteration, store the value and keep
        # iterating through `res`
        elif np.abs(target_avg_ani - current_avg_ani) <= np.abs(
            target_avg_ani - last_avg_weight
        ):
            last_res = res
            last_avg_weight = current_avg_ani
        # If the difference between the curent iteration average ANI and the target ANI
        # is greater than the last iteration, break the loop
        else:
            break
    return last_res


ani_graph = build_graph(snakemake.input.ANI, snakemake.input.FASTA)
if snakemake.params.leiden_resolution == "auto":
    leiden_resolution = pick_resolution(ani_graph)
else:
    leiden_resolution = snakemake.params.leiden_resolution
clusters = ani_graph.community_leiden(
    weights="weight", resolution_parameter=leiden_resolution
)

with open(snakemake.output[0], "w") as fout:
    fout.write(f"{leiden_resolution}\n")
    for i in reversed(np.argsort(clusters.sizes())):
        subgraph = clusters.subgraph(i)
        members = {v.attributes()["name"]: v.attributes()["length"] for v in subgraph.vs}
        members = sorted(members, key=members.get, reverse=True)
        fout.write(members[0] + "\t" + ",".join(members) + "\n")
