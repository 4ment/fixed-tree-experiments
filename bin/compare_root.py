import os
import re
import tempfile
from collections import defaultdict

import click
import dendropy


@click.command()
@click.argument("inputs", nargs=-1)
@click.option(
    "--outdir",
    type=click.Path(file_okay=False, dir_okay=True, path_type=str),
    default=".",
    help="output directory",
)
@click.option("--ref", type=click.File("r"), help="nexus file")
@click.option("--burnin", type=int, default=50001, help="burnin")
@click.option("--step", type=int, default=5, help="step size")
@click.option(
    "--chunk", type=int, default=500, help="number of trees to process at once"
)
def cli(inputs, outdir, ref, burnin, step, chunk):
    compare_root_txt = open(os.path.join(outdir, "compare_root.txt"), "w")
    dic = {}
    queries = {}
    taxa = dendropy.TaxonNamespace()

    for index, input_ in enumerate(inputs):
        counter = 0
        with open(input_, "r") as input1:
            for line in input1:
                line = line.strip()
                nexus = re.sub(r"\[[^\]]+\]", "", line)
                # replace index with name
                newick_list = list(filter(None, re.split(r"([,\(\):;])", nexus)))
                for idx, token in enumerate(newick_list):
                    if token not in ";:(,)" and newick_list[idx - 1] != ":":
                        if index == 0:
                            dic[token] = f"A{counter}A"
                            newick_list[idx] = dic[token]
                            counter += 1
                        else:
                            newick_list[idx] = dic[token]
                nexus = "".join(newick_list)

        tree = dendropy.Tree.get_from_string(
            nexus,
            schema="newick",
            rooting="default-rooted",
            taxon_namespace=taxa,
            preserve_underscores=True,
        )
        tree.encode_bipartitions()
        queries[input_] = frozenset(
            [e.split_bitmask for e in tree.seed_node.child_edge_iter()]
        )
        print(f"Parsed {input_} 1 tree with {len(taxa)} taxa")
        compare_root_txt.write(f"Parsed {input_} 1 tree with {len(taxa)} taxa\n")

        if index == 0:
            first_tree = tree

    tree_list = []
    header = []
    tree_counter = 0
    bitmasks = defaultdict(int)

    # with tempfile.NamedTemporaryFile(delete_on_close=False) as fw:
    for line in ref:
        line = line.strip()
        if line.endswith("Taxlabels"):
            # fw.write(line + "\n")
            header.append(line)
            for line in ref:
                line = line.strip()
                if line == ";":
                    # fw.write(line + "\n")
                    header.append(line)
                    break
                # fw.write("  " + dic[line.strip("'")] + "\n")
                header.append("  " + dic[line.strip("'").strip('"')])
        elif line.lower().endswith("translate"):
            # fw.write(line + "\n")
            header.append(line)
            for line in ref:
                line = line.strip()
                if line == ";":
                    # fw.write(line + "\n")
                    header.append(line)
                    break
                else:
                    line2 = line.rstrip(",")
                index, name = re.split(r"\s+", line2, maxsplit=1)
                name = name.strip("'").strip('"')
                end = "," if line.endswith(",") else ""
                # fw.write(index + " " + dic[name] + end + "\n")
                header.append(index + " " + dic[name] + end)
        elif line.startswith("tree"):
            tree_counter += 1
            if tree_counter < burnin or tree_counter % step != 0:
                continue
            tree_list.append(line)
            if len(tree_list) == chunk:
                with open("temp.nexus", "w") as fw:
                    fw.write("\n".join(header) + "\n")
                    fw.write("\n".join(tree_list) + "\nEnd;\n")
                trees = dendropy.TreeList.get_from_path(
                    "temp.nexus",
                    schema="nexus",
                    rooting="default-rooted",
                    taxon_namespace=taxa,
                    preserve_underscores=True,
                )
                for t in trees:
                    t.encode_bipartitions()
                    bitmasks[
                        frozenset(
                            [e.split_bitmask for e in t.seed_node.child_edge_iter()]
                        )
                    ] += 1

                tree_list.clear()
        elif len(tree_list) > 0:
            break
        else:
            header.append(line)
            # fw.write(line + "\n")

    if len(tree_list) != 0:
        with open("temp.nexus", "w") as fw:
            fw.write("\n".join(header) + "\n")
            fw.write("\n".join(tree_list) + "\nEnd;\n")
        trees = dendropy.TreeList.get_from_path(
            "temp.nexus",
            schema="nexus",
            rooting="default-rooted",
            taxon_namespace=taxa,
            preserve_underscores=True,
        )
        for t in trees:
            t.encode_bipartitions()
            bitmasks[
                frozenset([e.split_bitmask for e in t.seed_node.child_edge_iter()])
            ] += 1

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

    # Attempt to reroot the first tree using the best possible (highest posterior probability)
    # rooting from the reference trees. We only reroot if the posterior probability is greater than current
    # rooting.
    if inputs[0] in indices:
        found = False
        for i, tup in enumerate(sorted_bitmasks):
            if i == indices[inputs[0]]:
                break

            bitmask, count = tup
            for edge in first_tree.postorder_edge_iter():
                if edge.split_bitmask in bitmask:
                    first_tree.reroot_at_edge(edge, update_bipartitions=False)
                    reverse_dic = {v: k for k, v in dic.items()}
                    newick = first_tree.as_string(
                        schema="newick", suppress_internal_node_labels=True
                    )
                    newick = newick[newick.index("(") :].strip()
                    newick_list = list(filter(None, re.split(r"([,\(\):;])", newick)))
                    for idx, token in enumerate(newick_list):
                        if token not in ";:(,)" and newick_list[idx - 1] != ":":
                            newick_list[idx] = reverse_dic[token]
                    newick = "".join(newick_list)
                    with open(os.path.join(outdir, f"rerooted.tree"), "w") as fw:
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
