process combine_quast {
    publishDir "${params.final_outdir}", mode: 'copy', pattern: "*.yaml"
    tag "combine_quast_output"
    // debug true

    input:
    path(tsv_files)

    output:
    path ("combined_quast_mqc.yaml"), emit: yaml, optional:true


    shell:
    '''
    # comment here
    combine_quast.py -i "!{tsv_files}" \
                     -o combined_quast_mqc.yaml

    '''

}