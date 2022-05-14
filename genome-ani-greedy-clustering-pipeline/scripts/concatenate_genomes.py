from pathlib import Path


def is_binary(file_name):
    try:
        with open(file_name, 'tr') as check_file:
            check_file.read()
            return False
    except Exception:
        return True


input_path = Path(snakemake.input[0])
with open(snakemake.output[0], "w") as fout:
    for genome_path in input_path.iterdir():
        total_length = 0
        genome = genome_path.stem
        if not is_binary(genome_path):
            with open(genome_path) as fin:
                for line in fin:
                    if line.startswith(">"):
                        contig = f"{genome}_{total_length}"
                        fout.write(f">{contig}\n")
                    else:
                        fout.write(f"{line.strip()}\n")
                        total_length += len(line.strip())