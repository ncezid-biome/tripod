#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define the pipeline version
def pipeline_version = '1.0.0' // Replace with actual version or load dynamically

// Define the timestamp
def timestamp = new Date().format("yyyyMMdd_HHmmss")

// Define the output directory with the version and timestamp at runtime
params.final_outdir = params.outdir ? "${params.outdir}_v${pipeline_version}_${timestamp}" : "tripod_results_v${pipeline_version}_${timestamp}"
params.file_extension = "_v${pipeline_version}_${timestamp}"


Channel
    // search for pair-end wgs reads files in the given folder or any subfolders
   .fromFilePairs(["${params.wgs_reads}/*_R{1,2}*.fastq.gz", "${params.wgs_reads}/*_{1,2}*.fastq.gz", "${params.wgs_reads}/**/*_R{1,2}*.fastq.gz"], size: 2)
 .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
  .set { paired_wgs_reads }


Channel.fromPath(params.multiqc_config, checkIfExists: true).set { ch_config_for_multiqc }


include { shovill_assembly } from './modules/local/shovill_assembly.nf' 
include { quast } from './modules/quast/main.nf'
include { run_primersearch } from './modules/local/run_primersearch.nf' 
include { parse_primersearch } from './modules/local/parse_primersearch.nf' 
include { tripod_analysis as tripod_analysis_stool } from './modules/local/tripod_analysis.nf'
include { tripod_analysis as tripod_analysis_isolate } from './modules/local/tripod_analysis.nf'
include { combine_quast } from './modules/local/combine_quast.nf'
include { multiqc } from './modules/multiqc/main.nf' 
include { primer_performance } from './modules/local/primer_performance.nf' 

workflow {

    assembled_reads_ch = shovill_assembly(paired_wgs_reads)
    quast_report_ch = quast(assembled_reads_ch.fasta)
    assembled_reads_fasta = assembled_reads_ch.fasta.map { sample, fasta_path -> fasta_path }
    run_primersearch_ch = run_primersearch(assembled_reads_ch.fasta)
    extracted_amplicon_fasta = parse_primersearch(run_primersearch_ch.search_output)
    tripod_analysis_stool_ch = tripod_analysis_stool(extracted_amplicon_fasta.collect(), params.hmas_indir_stool, params.mapping_stool,'stool')
    tripod_analysis_isolate_ch = tripod_analysis_isolate(extracted_amplicon_fasta.collect(), params.hmas_indir_isolate, params.mapping_isolate,'isolate')
    combine_quast_ch = combine_quast(quast_report_ch.report.collect())
    primer_performance_ch = primer_performance(
      params.hmas_indir_stool,\
      params.hmas_indir_isolate,\
      params.good_sample_list_stool,\
      params.good_sample_list_isolate
    )

    multiqc(combine_quast_ch.yaml
            .combine(tripod_analysis_stool_ch.yaml)
            .combine(tripod_analysis_isolate_ch.yaml)
            .combine(primer_performance_ch.yaml), ch_config_for_multiqc)
}
