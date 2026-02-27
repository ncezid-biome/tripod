#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// -------------------------
// Pipeline metadata
// -------------------------
def pipeline_version = '1.0.0' 
def timestamp = new Date().format("yyyyMMdd_HHmmss")
params.final_outdir = params.outdir ? "${params.outdir}_v${pipeline_version}_${timestamp}" : "tripod_results_v${pipeline_version}_${timestamp}"
params.file_extension = "_v${pipeline_version}_${timestamp}"

// -------------------------
// Read mapping file and get allowed WGS IDs
// -------------------------
def allowed_ids = new File(params.mapping)
    .readLines()
    .drop(1) // skip header
    .collect { it.split(',') }
    .collect { it.size() > 2 ? it[2].trim() : "" } // convert missing to empty string
    .findAll { it } // remove empty strings
    as Set

def has_wgs = allowed_ids.size() > 0

// -------------------------
// Include modules
// -------------------------
include { shovill_assembly } from './modules/local/shovill_assembly.nf' 
include { run_primersearch } from './modules/local/run_primersearch.nf' 
include { parse_primersearch } from './modules/local/parse_primersearch.nf' 
include { tripod_analysis as tripod_analysis_stool } from './modules/local/tripod_analysis_oh_correlated.nf'
include { tripod_analysis as tripod_analysis_isolate } from './modules/local/tripod_analysis_oh_correlated.nf'
include { multiqc } from './modules/multiqc/main.nf' 
include { primer_performance } from './modules/local/primer_performance.nf' 
include { tripod_analysis as tripod_analysis_stool_isolate } from './modules/local/tripod_analysis_oh_correlated.nf'

// -------------------------
// MultiQC config channel
// -------------------------
Channel.fromPath(params.multiqc_config, checkIfExists: true).set { ch_config_for_multiqc }

// -------------------------
// WGS branch subworkflow
// -------------------------
workflow wgs_workflow {

    take:
    // nothing

    main:

    if (has_wgs) {

        Channel
            .fromFilePairs(
                [
                    "${params.wgs_reads}/*_R{1,2}*.fastq.gz",
                    "${params.wgs_reads}/*_{1,2}*.fastq.gz",
                    "${params.wgs_reads}/**/*_{1,2}*.fastq.gz",
                    "${params.wgs_reads}/**/*_R{1,2}*.fastq.gz"
                ],
                size: 2
            )
            .map { reads ->
                def sample_id = reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/, "")
                tuple(sample_id, reads[1])
            }
            .filter { sample_id, reads -> allowed_ids.contains(sample_id) }
            .set { paired_wgs_reads }

        assembled_reads_ch = shovill_assembly(paired_wgs_reads)
        run_primersearch_ch = run_primersearch(assembled_reads_ch.fasta)
        extracted_amplicon_fasta = parse_primersearch(run_primersearch_ch.search_output)

    } else {

        log.warn "No WGS IDs found — skipping assembly + primersearch."

        extracted_amplicon_fasta = Channel.empty()
    }

    emit:
    extracted_amplicon_fasta
}

// -------------------------
// Main workflow
// -------------------------
workflow {

    wgs_workflow()

    extracted_amplicon_fasta = wgs_workflow.out.extracted_amplicon_fasta

    extracted_amplicon_fasta
        .collect()
        .ifEmpty { [] }
        .set { collected_amplicons }

    tripod_analysis_stool_ch =
        tripod_analysis_stool(collected_amplicons, params.hmas_indir_stool, 'stool')

    tripod_analysis_isolate_ch =
        tripod_analysis_isolate(collected_amplicons, params.hmas_indir_isolate, 'isolate')

    tripod_analysis_stool_isolate_ch =
        tripod_analysis_stool_isolate(collected_amplicons, params.hmas_indir_stool, 'combined')

    def multiqc_input = tripod_analysis_stool_ch.yaml
                            .combine(tripod_analysis_isolate_ch.yaml)
                            .combine(tripod_analysis_stool_isolate_ch.yaml)

    multiqc(multiqc_input, ch_config_for_multiqc)
}