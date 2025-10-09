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
        # Handle partition format like "::3"
        slice_indices = slice(None, None, int(partition.split(":")[-1]))
    elif partition.startswith(":"):
        # Handle partition format like ":3:2" or ":3"
        slice_indices = slice(None, *list(map(int, partition.split(":")[1:])))
    elif "::" in partition:
        # Handle partition format like "3::2" or "3::"
        list_indices = partition.split("::")
        list_indices[0] = int(list_indices[0]) - 1
        list_indices[1] = None
        list_indices[2] = int(list_indices[2]) if list_indices[2] != "" else None
        slice_indices = slice(*list_indices)
    else:
        list_indices = list(map(int, partition.split(":")))
        list_indices[0] -= 1
        slice_indices = slice(*list_indices)

    with open(f"partition_{partition}.fasta", "w") as out_f:
        for header, sequence in sequences.items():
            out_f.write(f"{header}\n{sequence[slice_indices]}\n")
