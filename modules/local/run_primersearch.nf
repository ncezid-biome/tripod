process run_primersearch {
    publishDir "${params.final_outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    debug true

    input:
    tuple val(sample), path(fasta)


    output:
    tuple val(sample), path("${fasta.baseName}.ps"), path(fasta), emit: search_output, optional: true
    path "primersearch_error.log", emit: error, optional: true
    // path "versions.yml", emit: versions, optional: true

    shell:
    '''
    # Run primersearch
    if [[ -s "!{fasta}" ]] && [[ -s "!{params.primers}" ]]; then
        primersearch -seqall !{fasta} \
                    -infile !{params.primers} \
                    -mismatchpercent 6 \
                    -outfile !{fasta.baseName}.ps
    else
        echo "Error: Primersearch failed for !{fasta.baseName}" >> primersearch_error.log
    fi

    '''

}