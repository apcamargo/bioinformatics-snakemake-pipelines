import math
import os
import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from shutil import which


def count_sequences(file):
    n = 0
    with open(file) as fin:
        for line in fin:
            if line.startswith(">"):
                n += 1
    return n


def run_prodigal(workDir, currentId, chunkFile):
    cmd = [
        "prodigal",
        "-q",
        "-p",
        "meta",
        "-i",
        chunkFile.name,
        "-a",
        workDir.name + "/chunk" + str(currentId) + ".faa",
        "-o",
        workDir.name + "/chunk" + str(currentId) + ".out",
    ]
    subprocess.run(cmd, shell=False, check=True)


def append_fasta_file(file, startNum, targetFile):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    with open(targetFile, "a") as trgt:
        with open(file, "r") as input:
            for line in input:
                if line.startswith(">"):
                    match = re.match(pattern, line)
                    if match and match.group(3) == "1":
                        startNum = startNum + 1
                    line = (
                        match.group(1)
                        + str(startNum)
                        + "_"
                        + match.group(3)
                        + match.group(4)
                    )
                trgt.write(line.strip() + "\n")

if which("prodigal") is None:
    raise ValueError("Prodigal not found in the PATH.")

queryFile = snakemake.input[0]
proteinFile = snakemake.output[0]

tasks = max(snakemake.threads, 1)
n_sequences = count_sequences(queryFile)
seqsPerChunk = math.ceil(n_sequences / tasks)

# Delete opput file if it already exists. Otherwise the new results will be
# appended to it
if proteinFile and os.path.isfile(proteinFile):
    os.remove(proteinFile)

seqCnt = 0
currentChunk = 1

workDir = tempfile.TemporaryDirectory()
executor = ThreadPoolExecutor(max_workers=tasks)
currentFile = open(workDir.name + "/chunk" + str(currentChunk), "w")

with open(queryFile, "r") as fasta:
    for line in fasta:
        if line[0] == ">" and seqCnt == seqsPerChunk:
            currentFile.close()
            executor.submit(run_prodigal, workDir, currentChunk, currentFile)
            currentFile = None
            seqCnt = 0
            currentChunk += 1
        if currentFile is None:
            currentFile = open(workDir.name + "/chunk" + str(currentChunk), "w")
        currentFile.write(line)
        if line[0] == ">":
            seqCnt += 1
if seqCnt > 0:
    currentFile.close()
    executor.submit(run_prodigal, workDir, currentChunk, currentFile)

# await completion of tasks
executor.shutdown(wait=True)

# collect output
protIdStart = 0
for cur in range(1, currentChunk + 1):
    if proteinFile:
        append_fasta_file(
            workDir.name + "/chunk" + str(cur) + ".faa", protIdStart, proteinFile
        )
    protIdStart += seqsPerChunk
