#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.datasets = "${baseDir}/ssstree_datasets"
params.results = "results"


process RUN_IQTREE {
  publishDir "$params.results/datasets/${dataset}/iqtree", mode: 'copy'

  input: val(dataset)
  output:
    tuple val(dataset), path("ali.fasta.treefile"), env(SITES)
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

process CREATE_BEAST_XML{
  label 'ultrafast'

  publishDir "$params.results/datasets/${dataset}/${type}", mode: 'copy'

  input:
    tuple val(dataset), val(type), path(tree_file)
  output:
    tuple val(dataset), val(type), path("run01-fixed.xml")
  """
  helper.py beast --input ${params.datasets}/${dataset}/run01.xml --tree ${tree_file} --output run01-fixed.xml
  """
}

process RUN_BEAST{
  publishDir "$params.results/datasets/${dataset}/${type}", mode: 'copy'

  input:
    tuple val(dataset), val(type), path(xml_file)
  output:
    tuple val(dataset), val(type), path("*.log"), path("*.trees")
  """
  touch ${dataset}-${type}.txt
  $HOME/programs/BEASTv1.10.5pre_thorney_0.1.2/bin/beast -save_every 100000 -save_state checkpoint.state ${xml_file}
  """
}

process PARSE_MAP{
  label 'ultrafast'

  publishDir "$params.results/datasets/${dataset}/map", mode: 'copy'

  input:
    tuple val(dataset), path(log_file), path(trees_file)
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
    tuple val(dataset), path("${trees_file}-mcc.nexus"), path("${trees_file}-mcc.tree")
  """
  STATES=\$(grep tree ${trees_file}|wc -l)
  BURNIN=\$((STATES / 2))
  $HOME/programs/BEASTv1.10.5pre_thorney_0.1.2/bin/treeannotator -burninTrees \${BURNIN} ${trees_file} ${trees_file}-mcc.nexus
  python bin/tree2newick.py --input ${trees_file}-mcc.nexus --output ${trees_file}-mcc.tree
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
  
  RUN_IQTREE(datasets) // in: dataset out: ($dataset, $tree_file, $sites)
  RUN_LSD(RUN_IQTREE.out[0]) // in: $dataset, $tree_file, $sites out[0]: ($dataset, lsd2, $tree_file)
  RUN_TREETIME(RUN_IQTREE.out[0].map{it -> [it[0], it[1]]}) // in: $dataset, tree out[0]: ($dataset, treetime, $tree_file)
  
  ch_map = datasets.map{
    it ->
    [it, file("${params.datasets}/${it}/prob_run01.log"), file("${params.datasets}/${it}/run01.trees")]
  }
  PARSE_MAP(ch_map) // in: ($dataset, prob_run01.log, run01.trees) out: ($dataset, map, map.nwk)

  ch_mcc = datasets.map{
    it ->
    [it, "mcc", file("${params.datasets}/${it}/run01_burninNgen200000000_MCC.tree")]
  }

  ch_all = ch_mcc.concat(PARSE_MAP.out[0], RUN_LSD.out[0], RUN_TREETIME.out[0])
  CREATE_BEAST_XML(ch_all) // in: ($dataset, $type, $tree) out: ($dataset, $type, run01-fixed.xml)

  RUN_BEAST(CREATE_BEAST_XML.out) // in: ($dataset, $type, $xml) out: $dataset, $type, "*.log", "*.trees"
}
