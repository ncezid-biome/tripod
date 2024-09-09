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

// Channel
//     // search for pair-end HMAS reads files in the given folder or any subfolders
//    .fromFilePairs(["${params.hmas_reads}/*_R{1,2}*.fastq.gz", "${params.hmas_reads}/**/*_R{1,2}*.fastq.gz"], size: 2)
//  .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
//   .set { paired_hmas_reads }

Channel.fromPath(params.multiqc_config, checkIfExists: true).set { ch_config_for_multiqc }
// Channel.fromPath(params.custom_logo, checkIfExists: true).set { ch_logo_for_multiqc }


include { shovill_assembly } from './modules/local/shovill_assembly.nf' 
include { quast } from './modules/quast/main.nf'
include { extract_amplicon } from './modules/local/extract_amplicon.nf' 
include { tripod_analysis } from './modules/local/tripod_analysis.nf'
include { combine_quast } from './modules/local/combine_quast.nf'
include { multiqc } from './modules/multiqc/main.nf' 

workflow {
    // Filter out file pairs containing "Undetermined"
    // paired_reads = paired_reads.filter { pair -> 
    // !new File(pair[0]).getName().toLowerCase().startsWith("undetermined")}

    assembled_reads_ch = shovill_assembly(paired_wgs_reads)
    quast_report_ch = quast(assembled_reads_ch.fasta)
    assembled_reads_fasta = assembled_reads_ch.fasta.map { sample, fasta_path -> fasta_path }
    // quast_report_ch = quast(assembled_reads_fasta.collect())
    extracted_amplicon_ch = extract_amplicon(assembled_reads_ch.fasta)
    extracted_amplicon_fasta = extracted_amplicon_ch.fasta.map { sample, fasta_path -> fasta_path }
    tripod_analysis_ch = tripod_analysis(extracted_amplicon_fasta.collect())
    combine_quast_ch = combine_quast(quast_report_ch.report.collect())

    // // add fastqc and cutadapt log files (these are existing modules in MultiQC)
    // Channel.empty()
    //     .mix( FASTQC_RAW.out.fastqc_results )
    //     .mix( removed_primer_reads_ch.cutadapt_json )
    //     .map { sample, files -> files }
    //     .collect().ifEmpty([])
    //     .set { log_files }

    // // add custom content log files
    // multiqc(log_files
    //     .combine(ch_logo_for_multiqc)
    //     .combine(pear_log_ch.log)
    //     .combine(qfilter_log_ch.log)
    //     .combine(derep_log_ch)
    //     .combine(denoise_log_ch)
    //     .combine(combined_report_ch.primer_stats_mqc)
    //     .combine(combined_report_ch.read_length_mqc)
    //     .combine(combined_report_ch.report_mqc), ch_config_for_multiqc)

    // Channel.empty()
    //     .mix( quast_report_ch.report )
    //     .collect().ifEmpty([])
    //     .set { quast_out }

    multiqc(combine_quast_ch.yaml
            .combine(tripod_analysis_ch.yaml), ch_config_for_multiqc)
}
