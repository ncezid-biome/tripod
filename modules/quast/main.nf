process quast {
    // publishDir "${params.final_outdir}", mode: 'copy'
    tag "${sample}"

    input:
    tuple val(sample), path (fasta)

    output:
    path "report*.tsv"             , optional:true, emit: report

    shell:
    '''
    # comment here
    quast !{fasta} &&  mv quast_results/latest/report.tsv report_!{sample}.tsv

    '''

}