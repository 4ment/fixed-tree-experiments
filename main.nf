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
    tuple val(dataset), path("ali.fasta.lsd.date.nexus")
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

  publishDir "$params.results/datasets/${dataset}/fixed", mode: 'copy'

  input:
    tuple val(dataset), path(tree_file)
  output:
    tuple val(dataset), path("run1-fixed.xml")
  """
  helper.py beast --input ${params.datasets}/${dataset}/run1.xml --tree ${tree_file} --output run1-fixed.xml
  """
}

process RUN_BEAST{
  publishDir "$params.results/datasets/${dataset}/fixed", mode: 'copy'
  errorStrategy 'ignore'

  input:
    tuple val(dataset), path(xml_file)
  output:
    tuple val(dataset), path("*.log"), path("*.trees")
  """
  beast ${xml_file}
  """
}

def extract(path){
  m = path =~ /${params.datasets}\/(\w+)\/run1.xml/
  return m[0][1]
}

workflow {
  datasets = Channel.fromPath("${params.datasets}/**/run1.xml").map{
    it -> extract(it.toString())
  }
  
  RUN_IQTREE(datasets) // in: dataset out: (dataset, tree, sites)
  CREATE_DATE_FILE(datasets) //in: dataset out: (dataset, dates)
  RUN_LSD(RUN_IQTREE.out.join(CREATE_DATE_FILE.out)) // in: dataset, tree, sites, dates out[0]: (dataset, tree_file)
  CREATE_BEAST_XML(RUN_LSD.out[0]) // in: (dataset, tree_file) out: (dataset, xml)
  RUN_BEAST(CREATE_BEAST_XML.out) // in: (dataset, xml) out: dataset, "*.log", "*.trees"
}
