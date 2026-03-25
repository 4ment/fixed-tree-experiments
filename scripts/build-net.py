import os
import pickle
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import OrderedDict

import click
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from treezy import NewickReader, NexusReader, Tree
from treezy.tree_metric import RobinsonFouldsMetric


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


def canonicalize_value(bitset) -> int:
    return bitset.value if bitset[0] else (~bitset).value


def weighted_rf_matrix(
    beast_dicts: list[dict[int, float]],
    iqtree_dict: dict[int, float],
    mcc_dict: dict[int, float],
):
    beast_count = len(beast_dicts)
    n = len(beast_dicts) + 2
    mat = np.zeros((n, n))
    metric = RobinsonFouldsMetric()

    beast_rf = []
    iqtree_rf = []
    for i in range(n):
        for j in range(i + 1, n):
            # Determine which dict corresponds to index i
            if i < beast_count:
                dict_i = beast_dicts[i]
            elif i == beast_count:
                dict_i = iqtree_dict
            else:  # i == beast_count + 1
                dict_i = mcc_dict

            # Determine which dict corresponds to index j
            if j < beast_count:
                dict_j = beast_dicts[j]
            elif j == beast_count:
                dict_j = iqtree_dict
            else:
                dict_j = mcc_dict

            # Compute weighted RF
            mat[i, j] = mat[j, i] = metric.compute_from_splits(dict_i, dict_j)

            if i < beast_count and j < beast_count:
                beast_rf.append(mat[i, j])
            elif i < beast_count and j == beast_count:
                iqtree_rf.append(mat[i, j])

    beast_rf = np.array(beast_rf)
    iqtree_rf = np.array(iqtree_rf)
    z_score = (np.mean(iqtree_rf) - np.mean(beast_rf)) / np.std(beast_rf)
    print(beast_rf.shape, iqtree_rf.shape)
    print(z_score)

    max_dist = np.max(mat[np.triu_indices_from(mat, k=1)])
    similarity = 1 - mat / max_dist

    beast_rf = []
    iqtree_rf = []
    for i in range(n):
        for j in range(i + 1, n):
            if i < beast_count and j < beast_count:
                beast_rf.append(similarity[i, j])
            elif i < beast_count and j == beast_count:
                iqtree_rf.append(similarity[i, j])

    beast_rf = np.array(beast_rf)
    iqtree_rf = np.array(iqtree_rf)
    print((np.mean(iqtree_rf) - np.mean(beast_rf)) / np.std(beast_rf))

    return similarity, z_score


def rf_matrix(
    beast_sets: dict[set, int],
    iqtree_sets: dict[set, int],
    mcc_sets: dict[set, int],
    rb_sets: dict[set, int],
):
    beast_count = len(beast_sets)
    n = len(beast_sets) + len(rb_sets) + len(iqtree_sets) + len(mcc_sets)
    uniq_sets = (
        list(rb_sets.keys())
        + list(beast_sets.keys())
        + list(iqtree_sets.keys())
        + list(mcc_sets.keys())
    )

    mat = np.zeros((n, n))
    metric = RobinsonFouldsMetric()

    beast_rf = []
    iqtree_rf = []
    for i in range(n):
        for j in range(i + 1, n):
            mat[i, j] = mat[j, i] = metric.compute_from_splits(
                uniq_sets[i], uniq_sets[j]
            )
            if i < beast_count and j < beast_count:
                beast_rf.extend(
                    [mat[i, j]] * (beast_sets[uniq_sets[i]] * beast_sets[uniq_sets[j]])
                )
            elif i < beast_count and j == beast_count:
                iqtree_rf.extend(
                    [mat[i, j]] * (beast_sets[uniq_sets[i]] * iqtree_sets[uniq_sets[j]])
                )

    beast_rf = np.array(beast_rf)
    iqtree_rf = np.array(iqtree_rf)
    z_score = (np.mean(iqtree_rf) - np.mean(beast_rf)) / np.std(beast_rf)
    print(beast_rf.shape, iqtree_rf.shape)
    print(z_score)

    max_dist = np.max(mat[np.triu_indices_from(mat, k=1)])
    similarity = 1 - mat / max_dist

    beast_rf = []
    iqtree_rf = []
    for i in range(n):
        for j in range(i + 1, n):
            if i < beast_count and j < beast_count:
                beast_rf.extend(
                    [similarity[i, j]]
                    * (beast_sets[uniq_sets[i]] * beast_sets[uniq_sets[j]])
                )
            elif i < beast_count and j == beast_count:
                iqtree_rf.extend(
                    [similarity[i, j]]
                    * (beast_sets[uniq_sets[i]] * iqtree_sets[uniq_sets[j]])
                )

    beast_rf = np.array(beast_rf)
    iqtree_rf = np.array(iqtree_rf)
    print((np.mean(iqtree_rf) - np.mean(beast_rf)) / np.std(beast_rf))

    return similarity, z_score


def plot_graph(ax, G, colors, sizes, title=None):
    pos = nx.spring_layout(G, weight="weight", seed=42)
    node_colors = [colors[G.nodes[n]["method"]] for n in G.nodes]
    node_sizes = [sizes[G.nodes[n]["method"]] for n in G.nodes]
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, ax=ax)
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=10)


def build_graph(weighted, beast, iqtree, mcc, rb, rb_burnin, rb_step, dic):
    taxon_names = []
    rb_sets = defaultdict(int)
    beast_sets = defaultdict(int)
    iqtree_sets = defaultdict(int)
    mcc_sets = defaultdict(int)
    beast_dicts = []
    old2new = {}
    methods = {}
    rev_count = 0
    beast_count = 0
    iqtree_count = 0

    taxon_check = None

    if rb:
        if dic:
            for line in dic:
                new_name, name = line.strip().split(",")
                old2new[new_name] = name
            taxon_names = list(old2new.values())

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
                                frozenset(
                                    (
                                        n.descendant_bitset.value,
                                        (~n.descendant_bitset).value,
                                    )
                                )
                                # canonicalize_value(n.descendant_bitset)
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
                    beast_dict_blens = {
                        # canonicalize_value(n.descendant_bitset): n.distance
                        frozenset(
                            (n.descendant_bitset.value, (~n.descendant_bitset).value)
                        ): n.distance
                        for n in tree.nodes
                        if not n.is_root and not n.is_leaf
                    }
                    beast_dicts.append(beast_dict_blens)
                    beast_sets[frozenset(beast_dict_blens.keys())] += 1
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
        iqtree_dict_blens = {
            # canonicalize_value(n.descendant_bitset): n.distance
            frozenset(
                (n.descendant_bitset.value, (~n.descendant_bitset).value)
            ): n.distance
            for n in tree.nodes
            if not n.is_root and not n.is_leaf
        }
        bitsets = frozenset(iqtree_dict_blens.keys())
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
        mcc_dict_blens = {
            # canonicalize_value(n.descendant_bitset): n.distance
            frozenset(
                (n.descendant_bitset.value, (~n.descendant_bitset).value)
            ): n.distance
            for n in tree.nodes
            if not n.is_root and not n.is_leaf
        }
        mcc_sets[frozenset(mcc_dict_blens.keys())] += 1
        print(f"Number of MCC trees: {len(mcc_sets)} ({len(tree)})")
        methods.update(
            {
                i + rev_count + beast_count + iqtree_count: "MCC"
                for i in range(len(mcc_sets))
            }
        )

    if weighted:
        similarity, z_score = weighted_rf_matrix(
            beast_dicts, iqtree_dict_blens, mcc_dict_blens
        )
    else:
        similarity, z_score = rf_matrix(beast_sets, iqtree_sets, mcc_sets, rb_sets)

    # Convert to graph
    G = nx.Graph()

    # Add nodes
    n = similarity.shape[0]
    G.add_nodes_from(range(n))

    nx.set_node_attributes(G, methods, name="method")

    # Add weighted edges (optional: apply threshold)
    threshold = 0.1
    for i in range(n):
        for j in range(i + 1, n):
            # if similarity[i, j] > threshold:
            G.add_edge(i, j, weight=similarity[i, j])

    return G, z_score


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
    "--out",
    type=click.Path(dir_okay=False, writable=True, path_type=Path),
    help="Output file",
)
@click.option(
    "--multi",
    is_flag=True,
    help="Run all dataset",
)
@click.option(
    "--weighted",
    is_flag=True,
    help="Use Robinson-Foulds distance weighted by branch lengths",
)
@click.option(
    "--load", type=click.Path(exists=True, file_okay=False), help="Load from pickles"
)
@click.option(
    "--save",
    type=click.Path(exists=True, file_okay=False, writable=True),
    help="Save to pickles",
)
@click.option(
    "--results",
    type=click.Path(exists=False, file_okay=False, writable=True),
    default="results",
    help="Result folder generated by nextflow",
)
def cli(
    multi,
    weighted,
    beast,
    iqtree,
    mcc,
    rb,
    rb_burnin,
    rb_step,
    dic,
    out,
    load,
    save,
    results,
):
    if load and save:
        raise click.UsageError("Options --load and --save cannot be used together.")

    colors = {
        "BEAST": "#D3D3D3",
        # "RevBayes": "green",
        "IQ-TREE": "#009E73",
        "MCC": "#D55E00",
    }

    sizes = {"BEAST": 10, "MCC": 30, "IQ-TREE": 30}

    if not multi:
        G, z = build_graph(weighted, beast, iqtree, mcc, rb, rb_burnin, rb_step, dic)
        fig, ax = plt.subplots(figsize=(10, 10))
        plot_graph(ax, G, colors, sizes)
    else:
        datasets = [
            "fluH3L_bed15",
            "lassaL_kli22",
            "wnv_del20",
            "sars2_lem21",
            "ebov_dud17",
            "ebov_mba21",
            "fluH1L_bed15",
            "fluPb2_wor14",
            "fluVicL_bed15",
            "hiv_far14",
            "mumps_mon21",
            "rabies_via23",
            "sars2_can20",
            "sars2_pek22",
            "zikv_gru19",
        ]

        Gs = []
        zs = []

        for i, ds in enumerate(datasets):
            print(ds)

            if load is None:
                beast = open(f"{results}/datasets/{ds}/unconstrained/run01.trees", "r")
                iqtree = open(f"{results}/datasets/{ds}/iqtree/ali.fasta.treefile", "r")
                mcc = open(f"{results}/datasets/{ds}/mcc/mcc.tree", "r")
                dic_csv = None

                G, z = build_graph(weighted, beast, iqtree, mcc, None, -1, -1, dic_csv)
                Gs.append(G)
                zs.append(z)

                if save is not None:
                    with open(os.path.join(save, f"{ds}.pickle"), "wb") as f:
                        pickle.dump(G, f)
            else:
                with open(os.path.join(load, f"{ds}.pickle"), "rb") as f:
                    G = pickle.load(f)
                Gs.append(G)

        if save is not None:
            pd.DataFrame(zs, index=datasets, columns=["z_score"]).to_csv(
                os.path.join(save, "z_scores.csv")
            )
        elif load is not None:
            zs = pd.read_csv(os.path.join(load, "z_scores.csv"), index_col=0)[
                "z_score"
            ].tolist()

        pretty_labels = OrderedDict(
            {"BEAST": "BEAST Unconstrained", "MCC": "BEAST MCC", "IQ-TREE": "IQ-TREE"}
        )

        handles = [
            mpatches.Patch(color=colors[label], label=pretty_labels[label])
            for label in colors.keys()
        ]

        fig, axs = plt.subplots(4, 4, figsize=(15, 15))
        for i, idx in enumerate(np.argsort(zs)):
            print(datasets[idx], zs[idx])
            ax = axs[i // 4, i % 4]
            plot_graph(
                ax, Gs[idx], colors, sizes, title=f"{datasets[idx]}\nZ = {zs[idx]:.2f}"
            )

        axs[3, 3].axis("off")

        # place the legend in the empty subplot
        axs[3, 3].legend(
            handles=handles,
            loc="center",
            frameon=False,
            fontsize=12,
        )
        # fig.legend(handles=handles, title="Method", loc="center")
        plt.tight_layout()

    if out is not None:
        plt.savefig(out, bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    cli()
