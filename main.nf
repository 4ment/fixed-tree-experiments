#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.datasets = "${baseDir}/ssstree_datasets"
params.results = "results"

process RUN_COMPARE_ROOT {
  label 'ultrfast'

  publishDir "$params.results/datasets/${dataset}/", mode: 'copy'
  input:
    tuple val(dataset), path(treetime_tree), path(lsd_tree), path(ref)
  output:
    tuple dataset, path("compare_root.txt"), path("rerooted.tree", optional: true)

  """
  compare-root.py ${treetime_tree} ${lsd_tree} --ref ${ref}
  """
}

process RUN_REVBAYES {
  publishDir "$params.results/datasets/${dataset}/revbayes", mode: 'copy'

  input: val(dataset)
  output:
    val(dataset)
    path("revbayes.log")
    path("revbayes.trees")
    path("dic.csv")

  """
  helper.py rename --input ${params.datasets}/${dataset}/ali.fasta --output ali-renamed.fasta --dic_out dic.csv
  if [ "${dataset}" = "ebov_dud17" ]; then
    partition.py ali-renamed.fasta 1:14517 14518:18992
  fi
  rb ${params.datasets}/${dataset}/script.Rev
  """
}

process RUN_IQTREE {
  publishDir "$params.results/datasets/${dataset}/iqtree", mode: 'copy'

  input: val(dataset)
  output:
    tuple val(dataset), path("ali.fasta.treefile"), env(SITES)
    path("ali-renamed.fasta.iqtree")
    path("dic.csv")
    // env(RATE)
    // path("ali.fasta.timetree.nex")
    // path("ali.fasta.treefile")
    // path("ali.fasta.bionj")
    // path("ali.fasta.iqtree")
    // path("ali.fasta.log")
    // path("ali.fasta.mldist")
    // path("ali.fasta.uniqueseq.ph")
    // path("ali.fasta.timetree.lsd")
    // path("ali.fasta.timetree.lsd.nwk")// branch=subst

  """
  helper.py rename --input ${params.datasets}/${dataset}/ali.fasta --output ali-renamed.fasta --dic_out dic.csv
  iqtree2 -s ali-renamed.fasta -m GTR+G --prefix ali-renamed.fasta
  helper.py rename --input ali-renamed.fasta.treefile --output ali.fasta.treefile --dic_in dic.csv --reverse
  SITES=\$(grep "^Input data" ali-renamed.fasta.iqtree|awk '{print \$6}')
  """
}

process RUN_LSD {
  label 'ultrafast'

  publishDir "$params.results/datasets/${dataset}/lsd2", mode: 'copy'

  input:
    tuple val(dataset), path(tree_file),val(seq_length)
  output:
    tuple val(dataset), val("lsd2"), path("ali.fasta.lsd.date.nexus")
    env(RATE)
    path "ali.fasta.lsd.nwk", emit: lsd_tree_newick // branch=subst
    path "ali.fasta.lsd"
    path "dates.txt"

  """
  helper.py dates --input ${params.datasets}/${dataset}/run1.xml --output dates.txt
  lsd2 -i ${tree_file} \
       -d dates.txt \
       -o ali.fasta.lsd \
       -s ${seq_length} \
       -r a
  RATE=\$(grep "^ rate" ali.fasta.lsd|awk '{print \$2}'|sed "s/,//")
  """
}

process RUN_TREETIME {
  label 'ultrafast'

  publishDir "$params.results/datasets/${dataset}/treetime", mode: 'copy'

  input:
    tuple val(dataset),
          path(tree_file)
  output:
    tuple val(dataset), val("treetime"), path("timetree.tree")
    path "treetime/divergence_tree.nexus"
    path "treetime/molecular_clock.txt"
    path "treetime/sequence_evolution_model.txt"
    path "treetime/divergence_tree.nexus"
    path "treetime/timetree.nexus"
    path "dates.csv"

  """
  helper.py rename --input ${params.datasets}/${dataset}/ali.fasta --output ali-renamed.fasta --dic_out dic.csv
  helper.py rename --input ${tree_file} --output tree-renamed.tree --dic_in dic.csv
  helper.py dates --input ${params.datasets}/${dataset}/run1.xml --output dates.csv --dic dic.csv
  treetime --aln ali-renamed.fasta --tree tree-renamed.tree --dates dates.csv --outdir treetime
  helper.py rename --input treetime/timetree.nexus --output timetree.tree --dic_in dic.csv --reverse
  """
}

process RUN_BEAST_FIXED{
  publishDir "$params.results/datasets/${dataset}/${type}", mode: 'copy'

  input:
    tuple val(dataset), val(type), path(tree_file)
  output:
    tuple val(dataset), val(type), path("*.log"), path("*.trees")
    path("run01-fixed.xml")
  """
  helper.py beast --input ${params.datasets}/${dataset}/run01.xml --tree ${tree_file} --output run01-fixed.xml
  beast run01-fixed.xml
  """
}

process RUN_BEAST{
  publishDir "$params.results/datasets/${dataset}/unconstrained", mode: 'copy'

  input:
    val(dataset)
  output:
    tuple val(dataset), path("run01.trees")
    path("prob_run01.log")
    path("*.log")
  """
  sed 's/\\(chainLength="\\)[0-9]\\{1,\\}/\\11000/' ${params.datasets}/${dataset}/run01.xml > run01.xml
  beast run01.xml
  """
}

process PARSE_MAP{
  label 'ultrafast'

  publishDir "$params.results/datasets/${dataset}/map", mode: 'copy'

  input:
    tuple val(dataset), path(trees_file)
    path(log_file)
  output:
    tuple val(dataset), val("map"), path("map.nexus")
  """
  helper.py map --log ${log_file} --trees ${trees_file} > map.nexus
  """

}

process RUN_TREEANNOTATOR{
  publishDir "$params.results/datasets/${dataset}/mcc", mode: 'copy'

  input:
    tuple val(dataset), path(trees_file)
  output:
    tuple val(dataset), val("mcc"), path("mcc.tree")
    path("mcc.nexus")
  """
  treeannotator -burninTrees 200000000 ${trees_file} mcc.nexus
  helper.py nexus2newick --input mcc.nexus --output mcc.tree
  """
}

def extract(path){
  m = path =~ /${params.datasets}\/(\w+)\/run01.xml/
  return m[0][1]
}

workflow {
  datasets = Channel.fromPath("${params.datasets}/**/run01.xml").map{
    it -> extract(it.toString())
  }

  // Unconstrained BEAST run
  RUN_BEAST(datasets) // in: $dataset out: ($dataset, run01.trees), *log

  RUN_TREEANNOTATOR(RUN_BEAST.out[0]) // in: ($dataset, run01.trees) out: ($dataset, mcc, run01.trees-mcc.tree), run01.trees-mcc.nexus
  PARSE_MAP(RUN_BEAST.out[0], RUN_BEAST.out[1]) // in: ($dataset, run01.trees), prob_run01.log out: ($dataset, map, map.nwk)

  RUN_IQTREE(datasets) // in: dataset out: ($dataset, $tree_file, $sites)
  RUN_LSD(RUN_IQTREE.out[0]) // in: $dataset, $tree_file, $sites out[0]: ($dataset, lsd2, $tree_file)
  RUN_TREETIME(RUN_IQTREE.out[0].map{it -> [it[0], it[1]]}) // in: $dataset, tree out[0]: ($dataset, treetime, $tree_file)

  ch_compare_root = RUN_TREETIME.out[0].map{it -> [it[0], it[2]]}.join(RUN_LSD.out[0].map{it -> [it[0], it[2]]}).join(RUN_BEAST.out[0])
  RUN_COMPARE_ROOT(ch_compare_root) // in: ($dataset, $treetime_tree, $lsd_tree, $ref_trees) out: ($dataset, compare_root.txt, rerooted.tree)
  ch_rerooted = ch_compare_root.filter { dataset, compare_txt, maybe_file -> maybe_file.exists() }.map{it -> [it[0], "rerooted" it[2]]}
  ch_all = PARSE_MAP.out.join(RUN_LSD.out[0]).join(RUN_TREETIME.out[0]).join(RUN_TREEANNOTATOR.out[0]).join(ch_rerooted)
  RUN_BEAST_FIXED(ch_all) // in: ($dataset, $type, $tree) out: ($dataset, $type, run01-fixed.xml), *.log, *.trees
}
