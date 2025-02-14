process primer_performance {
    publishDir "${params.final_outdir}", mode: 'copy'
    tag "primer_performance"
    // debug true

    input:
    path(indir_stool)
    path(indir_isolate)
    path(good_sample_list_stool)
    path(good_sample_list_isoalte)

    output:
    path ("*.csv"), emit: report, optional:true
    path ("*.yaml"), emit: yaml, optional:true


    shell:
    '''
    get_primers_by_samples.py \
        -i !{indir_isolate} \
        -l !{indir_stool}\
        -is !{good_sample_list_isoalte}\
        -ls !{good_sample_list_stool}\

    '''

}