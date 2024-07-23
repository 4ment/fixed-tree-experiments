#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.datasets = "${baseDir}/ssstree_datasets"
params.results = "results"


process CREATE_DATE_FILE{
  label 'ultrafast'

  publishDir "$params.results/datasets/${dataset}/iqtree", mode: 'copy'

  input: val(dataset)
  output:
   tuple val(dataset), path("dates.txt")
  """
  helper.py dates --input ${params.datasets}/${dataset}/run1.xml --output dates.txt
  """
}

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
  iqtree2 -s ${params.datasets}/${dataset}/ali.fasta -m GTR+G --prefix ali.fasta
  SITES=\$(grep "^Input data" ali.fasta.iqtree|awk '{print \$6}')
  """
}

process RUN_LSD {
  label 'ultrafast'

  publishDir "$params.results/datasets/${dataset}/lsd2", mode: 'copy'

  input:
    tuple val(dataset),
          val(tree_file),
          val(seq_length),
          val(date_file)
  output:
    tuple val(dataset), val("lsd2"), path("ali.fasta.lsd.date.nexus")
    env(RATE)
    path "ali.fasta.lsd.nwk", emit: lsd_tree_newick // branch=subst
    path "ali.fasta.lsd"

  """
  lsd2 -i ${tree_file} \
       -d ${date_file} \
       -o ali.fasta.lsd \
       -s ${seq_length} \
       -r a
  RATE=\$(grep "^ rate" ali.fasta.lsd|awk '{print \$2}'|sed "s/,//")
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

def func(a){
  m = a =~ /results\/datasets\/(\w+)\/fixed\/run1-fixed.xml/
  return m[0][1]
}
def extract(path){
  m = path =~ /${params.datasets}\/(\w+)\/run01.xml/
  return m[0][1]
}

workflow {
  datasets = Channel.fromPath("${params.datasets}/**/run01.xml").map{
    it -> extract(it.toString())
  }
  // xml_files = files("${baseDir}/ssstree_datasets/**/run1.xml")
  // for(int i=0;i < xml_files.size; i++){
  //   m = xml_files[i] =~ /ssstree_datasets\/(\w+)\/run1.xml/
  //   xml_files[i] = m[0][1]
  // }
  
  RUN_IQTREE(datasets) // in: dataset out: (dataset, tree, sites)
  CREATE_DATE_FILE(datasets) //in: dataset out: (dataset, dates)
  RUN_LSD(RUN_IQTREE.out.join(CREATE_DATE_FILE.out)) // in: dataset, tree, sites, dates out[0]: (dataset, lsd2, $tree_file)
  all = Channel.from("fluH3L_bed15", "lassaL_kli22", "wnv_del20", "sars2_lem21", "ebov_dud17",
  "ebov_mba21", "fluH1L_bed15", "fluPb2_wor14", "fluVicL_bed15", "hiv_far14", "mumps_mon21", "rabies_via23", "sars2_can20",
  "sars2_pek22", "zikv_gru19")
  // all = Channel.from("ebov_mba21-1", "lassaL_kli22-1")
  // ch = all.map{
  //   it ->
  //   [it, file("${baseDir}/results/datasets/${it}/lsd2/ali.fasta.lsd.date.nexus")]
  // }

  ch_map = all.map{
    it ->
    [it, file("${params.datasets}/${it}/prob_run01.log"), file("${params.datasets}/${it}/run01.trees")]
  }
  PARSE_MAP(ch_map) // in: ($dataset, prob_run01.log, run01.trees) out: ($dataset, map, map.nwk)

  ch_mcc = all.map{
    it ->
    [it, "mcc", file("${params.datasets}/${it}/run01_burninNgen200000000_MCC.tree")]
  }

  ch_all = ch_mcc.concat(PARSE_MAP.out[0], RUN_LSD.out[0])
  CREATE_BEAST_XML(ch_all) // in: ($dataset, $type, $tree) out: ($dataset, $type, run01-fixed.xml)

  // xml_files = files("${baseDir}/results/datasets/**/fixed/run1-fixed.xml")
  // xml_ch = Channel.fromPath(xml_files)
  // ch = xml_ch.map{
  //   it ->
  //   [func(it.toString()), it]
  // }

  // ch_mcc = all.map{
  //   it ->
  //   [it, file("${params.datasets}/${it}/run1.trees")]
  // }
  // RUN_TREEANNOTATOR(ch_mcc) // in: (dataset, tree) out: (dataset, mcc_nexus, mcc_tree)

  // CREATE_BEAST_XML(ch)
  RUN_BEAST(CREATE_BEAST_XML.out) // in: (dataset, xml) out: dataset, "*.log", "*.trees"
}
