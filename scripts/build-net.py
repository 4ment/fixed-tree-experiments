from collections import defaultdict

import click
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from treezy import NewickReader, NexusReader, Tree
from treezy.tree_metric import RobinsonFouldsMetric
from pathlib import Path


def check_taxa(tree, taxon_check):
    if taxon_check is None:
        taxon_check = [None] * len(tree.taxon_names)
        for node in tree.nodes:
            if node.is_leaf:
                taxon_check[node.id] = node.name
    else:
        for node in tree.nodes:
            if node.is_leaf and node.name != taxon_check[node.id]:
                print(
                    f"Taxon {node.name} has different id {node.id} {taxon_check[node.id]}"
                )
                exit(1)
    return taxon_check


@click.command()
@click.option("--beast", type=click.File("r"), help="beast tree file")
@click.option("--iqtree", type=click.File("r"), help="iqtree tree file")
@click.option("--mcc", type=click.File("r"), help="MCC tree file")
@click.option(
    "--rb",
    type=click.Path(file_okay=True, dir_okay=False, path_type=str),
    help="revbayes tree file",
)
@click.option("--rb-burnin", type=int, default=1000, help="revbayes burnin")
@click.option("--rb-step", type=int, default=1, help="revbayes step size")
@click.option("--dic", type=click.File("r"), help="dic file")
@click.option(
    "--fig",
    type=click.Path(dir_okay=False, writable=True, path_type=Path),
    help="dic file",
)
def cli(beast, iqtree, mcc, rb, rb_burnin, rb_step, dic, fig):
    taxon_names = []
    rb_sets = defaultdict(int)
    beast_sets = defaultdict(int)
    iqtree_sets = defaultdict(int)
    mcc_sets = defaultdict(int)
    old2new = {}
    methods = {}
    rev_count = 0
    beast_count = 0
    iqtree_count = 0

    taxon_check = None

    if dic:
        for line in dic:
            new_name, name = line.strip().split(",")
            old2new[new_name] = name
        taxon_names = list(old2new.values())

    if rb:
        with open(rb, "r") as f:
            next(f)
            count = 0
            for idx, line in enumerate(f):
                if idx > rb_burnin and idx % rb_step == 0:
                    # if idx > 10000 and idx % 10 == 0:
                    line.strip()
                    newick = line[line.index("(") :]
                    tree = Tree.from_newick(newick)
                    if dic:
                        for node in tree.nodes:
                            if node.is_leaf:
                                node.name = old2new[node.name]
                        tree.taxon_names = taxon_names

                    taxon_check = check_taxa(tree, taxon_check)

                    # tree.make_binary()
                    tree.compute_descendant_bitset()
                    count += 1
                    rb_sets[
                        frozenset(
                            [
                                n.descendant_bitset.value
                                for n in tree.nodes
                                if not n.is_root and not n.is_leaf
                            ]
                        )
                    ] += 1
        rev_count = len(rb_sets)
        print(f"Number of revbayes trees: {count} # unique {rev_count} ({len(tree)})")
        methods.update({i: "RevBayes" for i in range(rev_count)})

    if beast:
        with NexusReader(beast, taxon_names, strip_quotes=True) as reader:
            counter = 0
            count = 0
            while reader.has_next():
                if counter > 10000 and counter % 100 == 0:
                    tree = reader.next()

                    taxon_check = check_taxa(tree, taxon_check)

                    tree.make_unrooted()
                    tree.compute_descendant_bitset()
                    beast_sets[
                        frozenset(
                            [
                                n.descendant_bitset.value
                                for n in tree.nodes
                                if not n.is_root and not n.is_leaf
                            ]
                        )
                    ] += 1
                    count += 1
                else:
                    reader.skip_next()
                counter += 1
        beast_count = len(beast_sets)
        print(f"Number of BEAST trees: {count} # unique {beast_count} ({len(tree)})")
        methods.update({i + rev_count: "BEAST" for i in range(beast_count)})

    if iqtree:
        with NewickReader(iqtree, taxon_names) as reader:
            tree = reader.next()

        taxon_check = check_taxa(tree, taxon_check)

        # tree.make_binary()
        tree.compute_descendant_bitset()
        bitsets = frozenset(
            [
                n.descendant_bitset.value
                for n in tree.nodes
                if not n.is_root and not n.is_leaf
            ]
        )
        iqtree_sets[bitsets] += 1
        iqtree_count = len(iqtree_sets)
        if bitsets in rb_sets:
            print(f"iqtree set in revbayes sets")
        else:
            print("iqtree set not in revbayes sets")
        print(f"Number of IQ-Tree trees: {1} # unique {iqtree_count} ({len(tree)})")
        methods.update(
            {i + rev_count + beast_count: "IQ-TREE" for i in range(len(iqtree_sets))}
        )

    if mcc is not None:
        with NexusReader(mcc, taxon_names, strip_quotes=True) as reader:
            tree = reader.next()

        taxon_check = check_taxa(tree, taxon_check)

        tree.make_unrooted()
        tree.compute_descendant_bitset()
        mcc_sets[
            frozenset(
                [
                    n.descendant_bitset.value
                    for n in tree.nodes
                    if not n.is_root and not n.is_leaf
                ]
            )
        ] += 1
        print(f"Number of MCC trees: {len(mcc_sets)} ({len(tree)})")
        methods.update(
            {
                i + rev_count + beast_count + iqtree_count: "MCC"
                for i in range(len(mcc_sets))
            }
        )

    n = len(beast_sets) + len(rb_sets) + len(iqtree_sets) + len(mcc_sets)
    uniq_sets = (
        list(rb_sets.keys())
        + list(beast_sets.keys())
        + list(iqtree_sets.keys())
        + list(mcc_sets.keys())
    )

    mat = np.zeros((n, n))
    metric = RobinsonFouldsMetric()
    for i in range(n):
        for j in range(i + 1, n):
            mat[i, j] = mat[j, i] = metric.compute_from_sets(uniq_sets[i], uniq_sets[j])

    max_dist = np.max(mat[np.triu_indices_from(mat, k=1)])
    similarity = 1 - mat / max_dist

    # Convert to graph
    G = nx.Graph()

    # Add nodes
    n = similarity.shape[0]
    G.add_nodes_from(range(n))

    nx.set_node_attributes(G, methods, name="method")
    colors = {
        "BEAST": "lightblue",
        "RevBayes": "green",
        "IQ-TREE": "orange",
        "MCC": "purple",
    }

    # Add weighted edges (optional: apply threshold)
    threshold = 0.1
    for i in range(n):
        for j in range(i + 1, n):
            # if similarity[i, j] > threshold:
            G.add_edge(i, j, weight=similarity[i, j])

    # Choose layout (2D)
    pos = nx.spring_layout(G, weight="weight", seed=42)
    # pos = nx.spectral_layout(G, weight="weight")

    # Draw network
    edges = G.edges(data=True)

    # weights = [d["weight"] for (_, _, d) in edges]
    # nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray', width=weights)
    # nx.draw_networkx_edge_labels(G, pos, edge_labels={(i,j): f"{d['weight']:.2f}" for i,j,d in edges})

    node_colors = [colors[G.nodes[n]["method"]] for n in G.nodes]

    # no edges
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=10)
    # nx.draw_networkx_labels(G, pos)

    handles = [
        mpatches.Patch(color=color, label=label) for label, color in colors.items()
    ]
    plt.legend(handles=handles, title="Method")

    if fig is not None:
        plt.axis("off")
        plt.savefig(fig, bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    cli()
