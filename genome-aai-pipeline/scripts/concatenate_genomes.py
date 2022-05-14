from pathlib import Path

input_path = Path(snakemake.input[0])
with open(snakemake.output[0], "w") as fout:
    for genome_path in input_path.iterdir():
        count = 1
        genome = genome_path.stem
        with open(genome_path) as fin:
            for line in fin:
                if line.startswith(">"):
                    contig = f"{genome}_{count}"
                    count += 1
                    fout.write(f">{contig}\n")
                else:
                    fout.write(f"{line.strip()}\n")