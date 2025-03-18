process tripod_analysis {
    publishDir "${params.final_outdir}", mode: 'copy'
    tag "tripod_analysis: ${tag_name}"
    // debug true

    input:
    path(fasta)
    path(indir)
    path(mapping)
    val(tag_name)

    output:
    path ("*.csv"), emit: report, optional:true
    path ("*.fasta"), emit: fasta, optional:true
    path ("*.yaml"), emit: yaml, optional:true


    shell:
    '''
    # comment here
    cat !{fasta} >> reference.fasta
    run_tripod_analysis.py \
        -i !{indir} \
        -r reference.fasta \
        -p !{mapping} \
        -o "tripod_!{tag_name}_!{params.file_extension}.csv" \
        -t !{tag_name}

    '''

}