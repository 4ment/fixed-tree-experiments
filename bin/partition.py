#!/usr/bin/env python

import sys

partitions = [arg for arg in sys.argv[2:]]
print(f"Partitions: {partitions}")

header = None
sequence = None
sequences = OrderedDict = {}

with open(sys.argv[1], "r") as f:
    for line in f:
        if line.startswith(">"):
            if sequence:
                sequences[header] = sequence
            header = line.strip()
            sequence = ""
        else:
            sequence += line.strip()
if sequence:
    sequences[header] = sequence

for partition in partitions:
    if partition.startswith("::"):
        slice_indices = slice(None, None, int(partition.split(":")[-1]))
    elif partition.startswith(":"):
        slice_indices = slice(None, *list(map(int, partition.split(":")[1:])))
    else:
        list_indices = list(map(int, partition.split(":")))
        list_indices[0] -= 1
        slice_indices = slice(*list_indices)

    with open(f"partition_{partition}.fasta", "w") as out_f:
        for header, sequence in sequences.items():
            out_f.write(f"{header}\n{sequence[slice_indices]}\n")
