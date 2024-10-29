#!/usr/bin/env python

import re
import xml.etree.ElementTree as ET

import click


def convert(fp):
    dic = {}
    for line in fp:
        line = line.strip()
        if line.lower().endswith("translate"):
            for line in fp:
                line = line.strip()
                if line == ";":
                    break
                else:
                    line2 = line.rstrip(",")
                index, name = re.split(r"\s+", line2, maxsplit=1)
                # dic[int(index)] = name
                dic[index] = name
        elif line.startswith("tree"):
            start = line.index("(")
            nexus = line[start:]
            # remove comments such as [%height=1]
            nexus = re.sub(r"\[[^\]]+\]", "", nexus)
            # replace index with name
            if len(dic) > 0:
                newick_list = list(filter(None, re.split(r"([,\(\):;])", nexus)))
                for idx, token in enumerate(newick_list):
                    if token not in ";:(,)" and newick_list[idx - 1] != ":":
                        newick_list[idx] = dic[token]
                nexus = "".join(newick_list)
                # for index in sorted(dic.keys(), reverse=True):
                #     name_escaped = re.escape(dic[index])
                #     print(rf"([,\(]){index}:", r"\1" + name_escaped + ":")
                #     nexus = re.sub(r"([,\(])" + f"{index}:", r"\1" +dic[index] + ":", nexus, count=1)
                #     print(nexus)
                #     # nexus = nexus.replace(f"{index}:", f"{dic[index]}:")
            return nexus


def renamer(id_):
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
@click.option("--rename", default=False)
def beast(input, output, tree, rename):
    with open(tree, "r") as fp:
        for line in fp:
            if line.startswith("#NEXUS"):
                newick = convert(fp)
            else:
                newick = line.strip()
            break

    tree = ET.parse(input)
    root = tree.getroot()

    # divide chain length by 10
    mcmc_elem = root.find("mcmc")
    mcmc_elem.attrib[
        "chainLength"
    ] = f"{int(int(mcmc_elem.attrib['chainLength']) / 10)}"

    logs_elem = mcmc_elem.findall("log")
    logs_to_remove = ("sitellLog", "opsTopoLog")
    for log_elem in logs_elem:
        if log_elem.attrib["id"] in logs_to_remove:
            mcmc_elem.remove(log_elem)

    logtrees_elem = mcmc_elem.findall("logTree")
    for logtree_elem in logtrees_elem:
        if logtree_elem.attrib["id"] == "topoFileLog":
            mcmc_elem.remove(logtree_elem)

    logCladeOperated = mcmc_elem.find("logCladeOperated")
    if logCladeOperated is not None:
        mcmc_elem.remove(logCladeOperated)

    root.remove(root.find("siteLogLikelihood"))

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

    if rename:
        all_taxa = {}
        for taxa in root.findall("taxa"):
            for taxon in taxa:
                if taxon.get("id") is not None:
                    id_ = taxon.get("id")
                    all_taxa[taxon.get("id")] = renamer(id_)

        for taxon, renamed in all_taxa.items():
            if taxon != renamed:
                if " " in taxon:
                    taxon = '"' + taxon + '"'
                newick = newick.replace(renamed, taxon)

    newick_element = ET.Element("newick", id="startingTree")
    newick_element.text = newick
    root.insert(index, newick_element)

    tree.write(output)


@click.command(help="Create date file for LSD or treetime")
@click.option("--input", type=click.UNPROCESSED, required=True, help="BEAST XML file")
@click.option(
    "--output", type=click.UNPROCESSED, required=True, help="file containing dates"
)
@click.option(
    "--dic", type=click.UNPROCESSED, required=False, help="file renaming sequences"
)
def dates(input, output, dic):
    tree = ET.parse(input)
    root = tree.getroot()
    taxa = root.find("taxa")
    old2new = None
    if dic is not None:
        old2new = {}
        with open(dic, "r") as fp:
            for line in fp:
                new_name, name = line.strip().split(",")
                old2new[name] = new_name
    # else:
    #     for taxon in taxa:
    #         names.append(taxon.get("id"))
    #         dates[taxon.get("id")] = taxon.find("date").get("value")
    #     counter = 0
    #     for name in sorted(names):
    #         old2new[name] = f"A{counter}A"
    #         counter += 1

    with open(output, "w") as fp:
        if output.endswith("csv"):
            fp.write("name, date\n")
        else:
            fp.write(str(len(taxa)) + "\n")
        for taxon in taxa:
            id_ = taxon.get("id")
            if old2new is not None:
                id_ = old2new[id_]

            # id_ = renamer(id_)
            date = taxon.find("date").get("value")
            if output.endswith("csv"):
                fp.write(f"{id_},{date}\n")
            else:
                fp.write(f"{id_}\t{date}\n")


@click.command(help="Exctract MAP tree in BEAST output")
@click.option("--log", type=click.File("r"), help="log file containing posterior trace")
@click.option("--trees", type=click.File("r"), help="log file containing trees")
@click.option("--output", type=click.File("w"), default=click.get_text_stream("stdout"))
def map(log, trees, output):
    max_prob = -float("inf")
    for line in log:
        if line.startswith("#"):
            continue
        elif line.startswith("state"):
            idx_posterior = line.split("\t").index("posterior")
        else:
            a = line.split("\t")
            prob = float(float(a[idx_posterior].strip()))
            s = a[0]
            if int(s) % 10000 == 0 and prob > max_prob:
                max_prob = prob
                state = s
    # print(state, max_prob)

    for line in trees:
        if line.startswith("tree"):
            a = re.split(r"\s+", line)
            s = a[1].split("_")[1]
            if s == state:
                output.write(line)
        else:
            output.write(line)


@click.command(help="Rename fasta or tree files")
@click.option("--input", type=click.UNPROCESSED, required=True, help="input file")
@click.option(
    "--output",
    type=click.UNPROCESSED,
    required=True,
    help="output file",
)
@click.option("--dic_in", type=click.File("r"), help="file for dictionary")
@click.option("--dic_out", type=click.File("w"), help="file for dictionary")
@click.option("--reverse", is_flag=True)
def rename(input, output, dic_in, dic_out, reverse):
    dictionary = None
    if dic_in:
        dictionary = {}
        for line in dic_in:
            if "," in line:
                new_name, name = line.strip().split(",")
                if reverse:
                    dictionary[new_name] = name
                else:
                    dictionary[name] = new_name
    # print(dictionary)

    is_seq = False
    with open(input, "r") as fp:
        for line in fp:
            if line.startswith(">"):
                is_seq = True
            break

    if is_seq:
        # rename fasta file sequences
        counter = 0
        names = []
        dic = {}
        with open(input, "r") as fp:
            for line in fp:
                if line.startswith(">"):
                    names.append(line[1:].strip())

        for name in sorted(names):
            dic[name] = f"A{counter}A"
            dic_out.write(f"A{counter}A,{name}\n")
            counter += 1

        with open(output, "w") as out:
            with open(input, "r") as fp:
                for line in fp:
                    if line.startswith(">"):
                        header = line[1:].strip()
                        out.write(f">{dic[header]}\n")
                    else:
                        out.write(line)
    else:
        with open(input, "r") as fp:
            for line in fp:
                line = line.strip()
                # Tree tree1 is for treetime (nexus file)
                if line.startswith("(") or line.startswith("Tree tree1"):
                    if line.startswith("Tree tree1"):
                        start = line.index("(")
                        line = line[start:]
                    # remove comments such as [&height=1]
                    line = re.sub(r"\[[^\]]+\]", "", line)
                    newick_list = list(filter(None, re.split(r"([,\(\):;])", line)))
                    for idx, token in enumerate(newick_list):
                        if (
                            token not in ";:(,)"
                            and newick_list[idx - 1] != ":"
                            and newick_list[idx - 1] != ")"
                        ):
                            newick_list[idx] = dictionary[token]

        with open(output, "w") as fp:
            fp.write("".join(newick_list))


@click.group(help="CLI tool to convert and generate files")
def cli():
    pass


cli.add_command(beast)
cli.add_command(dates)
cli.add_command(map)
cli.add_command(rename)

if __name__ == "__main__":
    cli()
