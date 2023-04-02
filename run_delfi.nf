/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DELFI NF pipeline : 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : 
    Author : Apiwat Sangphukieo
    Email  : sangphukieo@gmail.com
    ref    :
    usage  : nextflow run run_delfi.nf -with-report -with-dag dag.png -with-docker asangphukieo/xdelfi:v1.0.3 -c test.config -profile docker_test
----------------------------------------------------------------------------------------
*/

// Declare syntax version
nextflow.enable.dsl=2
// Script parameters
params.query = "/home/apiwat/hdd/DELFI/delfi_scripts/INPUT/sample/*.bgz"
params.lib = "/home/apiwat/hdd/DELFI/delfi_scripts/DELFI/iDELFI.R"
params.exome = "Agilent_Human_Exon_V6_UTRs"
params.hg = "hg19"
params.metadata = "/home/apiwat/hdd/DELFI/delfi_scripts/DELFI_NF/metadata_all_2235.csv"
params.script = "/home/apiwat/hdd/DELFI/delfi_scripts/DELFI_NF/script"
params.publish_dir = "/home/apiwat/hdd/DELFI/delfi_scripts/DELFI_NF/TEST2"

log.info """\
 D E L F I - N F   P I P E L I N E
 ===================================
 query            : ${params.query}
 libraly          : ${params.lib}
 exome prob       : ${params.exome}
 ref version      : ${params.hg}
 metadata file    : ${params.metadata}
 script path      : ${params.script}
 output folder    : ${params.publish_dir}
 """

process read_fragmentGC {
  input:
    val script
    path query
    path lib
    val exome
  output:
    path "./02_fragmentGC/${query.getSimpleName()}_frags.rds", emit: rds

    """
    Rscript $script/read_fragmentGC.r -f $query -l $lib -o ./02_fragmentGC -e $exome
    """
}

process ab_bin {
  input:
    val script
    val hg
    path lib
  output:
    path "AB.rds", emit: rds
    path "filters.hg19.rda", emit: rda
    
    """
    Rscript $script/ab_bin.r -g $hg -l $lib
    """
}

process filter_fragment_compartment_100kb {
  input:
    val script
    path tsv_frags
    path lib
    val ab
    val hg
    path filter
  output:
    path "./03_count_100kb/${tsv_frags.getSimpleName()}_count_100kb.rda", emit: count_100kb

    """
    Rscript $script/filter_fragment_compartment_100kb.r -f $tsv_frags -l $lib -g $hg -a $ab -b $filter
    """

}

process gc_correct_compartment_100kb {
  input:
    val script
    path count_100kb
    path lib
  output:
    path "./04_bins_100kb/${count_100kb.getSimpleName()}_bin_100kb.rds", emit: bin_100kb

    """
    Rscript $script/gc_correct_compartment_100kb.r -f $count_100kb -l $lib
    """  

}

process combine_bin_compartment {
  input:
    val script
    path bin_100kb
    path lib
  output:
    path "05_combine_100kb_bin/bins_100kbcompartments.rds", emit: bins_100kbcompartments

    """
    Rscript $script/combine_bin_compartment.r -d . -l $lib
    """  

}

process bin_5mb {
  publishDir "${params.publish_dir}", mode: 'copy', overwrite: true
  input:
    val script
    path bins_100kbcompartments
    path lib
  output:
    path "06_combine_5mb_bin/bins_5mbcompartments.rds", emit: bins_5mbcompartments

    """
    Rscript $script/bin_5mb.r -f $bins_100kbcompartments -l $lib
    """  
}

process bin_chr_arm {
  publishDir "${params.publish_dir}", mode: 'copy', overwrite: true
  input:
    val script
    path bins_100kbcompartments
    path lib
  output:
    path "06_02_combine_chr/bins_chr_compartments.rds", emit: bins_chr

    """
    Rscript $script/bin_chr_arm.r -f $bins_100kbcompartments -l $lib
    """  
}

process meta_data {
  publishDir "${params.publish_dir}", mode: 'copy', overwrite: true
  input:
    val script
    path bins_5mbcompartments
    path lib
    path metadata
  output:
    path "07_summarize_data/summary_data.csv.rds", emit: sum_data

    """
    Rscript $script/meta_data.r -f $bins_5mbcompartments -l $lib -t $metadata
    """  
}

process build_feature {
  publishDir "${params.publish_dir}", mode: 'copy', overwrite: true
  input:
    val script
    path bins_5mbcompartments
    path bins_chr
    path lib
    path sumdata
  output:
    path "*"

    """
    Rscript $script/build_feature.r -f $bins_5mbcompartments -b $bins_chr -l $lib -t $sumdata
    """  
}

workflow {
   def query_ch = Channel.fromPath(params.query)
   def script_ch = Channel.fromPath(params.script)
   def lib_ch = Channel.fromPath(params.lib) 
   def metadata_ch = Channel.fromPath(params.metadata)  
   read_fragmentGC(script_ch,query_ch, lib_ch, params.exome) 
   ab_bin(script_ch,params.hg,lib_ch)
   filter_fragment_compartment_100kb(script_ch, read_fragmentGC.out.rds, lib_ch, ab_bin.out.rds,params.hg, ab_bin.out.rda)
   gc_correct_compartment_100kb(script_ch, filter_fragment_compartment_100kb.out.count_100kb,lib_ch)
   combine_bin_compartment(script_ch, gc_correct_compartment_100kb.out.collect(),lib_ch)
   bin_5mb(script_ch, combine_bin_compartment.out.bins_100kbcompartments,lib_ch)
   meta_data(script_ch,bin_5mb.out.bins_5mbcompartments,lib_ch,metadata_ch)
   bin_chr_arm(script_ch,combine_bin_compartment.out.bins_100kbcompartments,lib_ch)
   build_feature(script_ch,bin_5mb.out.bins_5mbcompartments,bin_chr_arm.out.bins_chr,lib_ch,meta_data.out.sum_data )

}

