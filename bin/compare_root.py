#!/usr/bin/env python

import os
from collections import defaultdict

import click
from treezy import NewickReader, NexusReader


@click.command(
    help="""
Compare rootings of trees against a reference tree file.

This script reads multiple input tree files and a reference tree file in Nexus format,
and compares the rootings of the input trees against the rootings in the reference
tree. It counts how many times each unique rooting appears in the reference tree,
and checks if the first input tree can be rerooted to match the best rooting from the
reference tree.
"""
)
@click.argument("inputs", nargs=-1)
@click.option(
    "--outdir",
    type=click.Path(file_okay=False, dir_okay=True, path_type=str),
    default=".",
    help="output directory",
)
@click.option("--ref", type=click.File("r"), help="reference trees in nexus file")
@click.option("--burnin", type=int, default=50001, help="burnin")
@click.option("--step", type=int, default=5, help="step size")
def cli(inputs, outdir, ref, burnin, step):
    compare_root_txt = open(os.path.join(outdir, "compare_root.txt"), "w")
    queries = {}
    taxa = []

    for input_ in inputs:
        with NewickReader(input_, taxa) as reader:
            tree = reader.next()
        tree.compute_descendant_bitset()

        queries[input_] = frozenset(
            {c.descendant_bitset.value for c in tree.root.children}
        )

        print(f"Parsed {input_} 1 tree with {len(taxa)} taxa")
        compare_root_txt.write(f"Parsed {input_} 1 tree with {len(taxa)} taxa\n")

    tree_counter = 0
    bitmasks = defaultdict(int)

    with NexusReader(ref, taxa) as reader:
        while reader.has_next():
            tree_counter += 1

            if tree_counter >= burnin and tree_counter % step == 0:
                tree = reader.next()
                tree.compute_descendant_bitset()
                bitmasks[
                    frozenset({c.descendant_bitset.value for c in tree.root.children})
                ] += 1
            else:
                reader.skip_next()

    tree_count = sum(bitmasks.values())

    print(
        f"Parsed {ref.name} {tree_count} trees ({tree_counter}) with {len(taxa)} taxa\n"
    )
    compare_root_txt.write(
        f"Parsed {ref.name} {tree_count} trees ({tree_counter}) with {len(taxa)} taxa\n\n"
    )

    for filename, b in queries.items():
        count = bitmasks.get(b, 0)
        print(f"{filename} {count} {tree_count} {len(bitmasks)}")
        compare_root_txt.write(f"{filename} {count} {tree_count} {len(bitmasks)}\n")

    print()
    compare_root_txt.write("\n")

    sorted_bitmasks = sorted(bitmasks.items(), key=lambda x: x[1], reverse=True)
    indices = {}
    for i, tup in enumerate(sorted_bitmasks):
        bitmask, count = tup
        m = ""
        for filename, b in queries.items():
            if b == bitmask:
                indices[filename] = i
                m += f"{inputs.index(filename)+1} "

        print(f"{i+1} {count} {count/tree_count} {m}")
        compare_root_txt.write(f"{i+1} {count} {count/tree_count} {m}\n")

    # Attempt to reroot the first tree using the best possible
    # (highest posterior probability) rooting from the reference trees.
    # We only reroot if the posterior probability is greater than current rooting.
    if inputs[0] in indices:
        with NewickReader(inputs[0], taxa) as reader:
            tree = reader.next()
        tree.compute_descendant_bitset()

        found = False
        for i, tup in enumerate(sorted_bitmasks):
            if i == indices[inputs[0]]:
                break

            bitmask, count = tup
            for node in tree.postorder():
                if node.descendant_bitset.value in bitmask:
                    tree.reroot_above(node)
                    newick = tree.newick()
                    with open(os.path.join(outdir, "rerooted.tree"), "w") as fw:
                        fw.write(newick + "\n")

                    print(f"\nRerooted from {indices[inputs[0]]+1} to {i+1}")
                    compare_root_txt.write(
                        f"\nRerooted from {indices[inputs[0]]+1} to {i+1}\n"
                    )
                    found = True
                    break
            if found:
                break

    compare_root_txt.close()


if __name__ == "__main__":
    cli()
