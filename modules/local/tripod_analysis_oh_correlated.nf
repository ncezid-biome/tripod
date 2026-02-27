process tripod_analysis {
    publishDir "${params.final_outdir}", mode: 'copy'
    tag "tripod_analysis: ${tag_name}"
    // debug true

    input:
    path(fasta)
    path(indir)
    val(tag_name)

    output:
    path ("*.csv"), emit: report, optional:true
    path ("*.fasta"), emit: fasta, optional:true
    path ("*.yaml"), emit: yaml, optional:true


    shell:
    '''
    # comment here
    cat !{fasta} >> reference.fasta

    # Conditionally run the 2nd script if tag_name == "combined"
    if [ "!{tag_name}" = "combined" ]; then
        run_tripod_combined_analysis_oh_correlated.py \
            -i !{indir} \
            -r reference.fasta \
            -s !{params.hmas_indir_isolate} \
            -o "combined_output_!{params.file_extension}.csv" \
            -t !{tag_name} \
            -c !{params.mapping}
    else
        run_tripod_analysis_oh_correlated.py \
            -i !{indir} \
            -r reference.fasta \
            -o "tripod_!{tag_name}_!{params.file_extension}.csv" \
            -t !{tag_name} \
            -c !{params.mapping}
    fi

    '''

}