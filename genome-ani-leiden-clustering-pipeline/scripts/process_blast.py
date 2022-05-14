with open(snakemake.output[0], "w") as fout:
    with open(snakemake.input[0]) as fin:
        for line in fin:
            line = line.split()
            qoffset = int(line[0].rsplit("_", 1)[1])
            toffset = int(line[1].rsplit("_", 1)[1])
            line[0] = line[0].rsplit("_", 1)[0]
            line[1] = line[1].rsplit("_", 1)[0]
            line[4] = str(int(line[4]) + qoffset)
            line[5] = str(int(line[5]) + qoffset)
            line[6] = str(int(line[6]) + toffset)
            line[7] = str(int(line[7]) + toffset)
            fout.write("\t".join(line) + "\n")