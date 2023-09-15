#!/usr/bin/env python

import re
import xml.etree.ElementTree as ET

import click


def rename(id_):
    # rename sequences because iqtree changes some symbols
    # zikv_gru19: one sequence has \ in its name
    # rabies_via23: some sequence names contain @
    # fluPb2_wor14: *
    # fluH3L_bed15: +
    # wnv_del20: '
    # if there is space inside the name, iqtree truncates the string up to the space
    if " " in id_:
        id_ = id_[0 : id_.index(" ")]
    return re.sub(r"[@\\\*\?\+']", "_", id_)


@click.command(help="Create BEAST file with fixed topology")
@click.option("--input", type=click.UNPROCESSED, required=True, help="BEAST XML file")
@click.option(
    "--output",
    type=click.UNPROCESSED,
    required=True,
    help="BEAST XML file with fixed topology",
)
@click.option("--tree", type=click.UNPROCESSED, required=True, help="tree file")
def beast(input, output, tree):
    with open(tree, "r") as fp:
        for line in fp:
            if line.startswith("tree 1 = "):
                newick = line.replace("tree 1 = ", "").strip()
    newick = re.sub(r'\[&date="\d+\.?\d*"]', "", newick)

    tree = ET.parse(input)
    root = tree.getroot()

    # remove operators that modify the topology
    operators = ("subtreeSlide", "narrowExchange", "wideExchange", "wilsonBalding")

    operators_elem = root.find("operators")
    for operator_str in operators:
        for operator in operators_elem.findall(operator_str):
            operators_elem.remove(operator)

    # fix topology
    simulator = root.find("coalescentSimulator")
    if simulator is None:
        simulator = root.find("coalescentTree")

    index = list(root).index(simulator)
    root.remove(simulator)

    all_taxa = {}
    for taxa in root.findall("taxa"):
        for taxon in taxa:
            if taxon.get("id") is not None:
                id_ = taxon.get("id")
                all_taxa[taxon.get("id")] = rename(id_)

    for taxon, renamed in all_taxa.items():
        if taxon != renamed:
            if " " in taxon:
                taxon = '"' + taxon + '"'
            newick = newick.replace(renamed, taxon)
    newick_element = ET.Element("newick", id="startingTree")
    newick_element.text = newick
    root.insert(index, newick_element)

    tree.write(output)


@click.command(help="Create date file for LSD")
@click.option("--input", type=click.UNPROCESSED, required=True, help="BEAST XML file")
@click.option(
    "--output", type=click.UNPROCESSED, required=True, help="file containing dates"
)
def dates(input, output):
    tree = ET.parse(input)
    root = tree.getroot()
    taxa = root.find("taxa")
    with open(output, "w") as fp:
        fp.write(str(len(taxa)) + "\n")
        for taxon in taxa:
            id_ = taxon.get("id")
            id_ = rename(id_)
            date = taxon.find("date").get("value")
            fp.write(f"{id_}\t{date}\n")


@click.group(help="CLI tool to convert and generate files")
def cli():
    pass


cli.add_command(beast)
cli.add_command(dates)

if __name__ == "__main__":
    cli()
