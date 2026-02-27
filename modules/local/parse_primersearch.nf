process parse_primersearch {
    publishDir "${params.final_outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    debug true

    input:
    tuple val(sample), path(ps_file), path(fasta)

    output:
    path ("**/*_extractedAmplicons.fasta"), emit: fasta, optional:true

    shell:
    '''
    parse_primersearch.py --sequence !{fasta} \
                          --primersearch !{ps_file} \
                          

    #cat <<-END_VERSIONS > versions.yml
    #"${task.process}":
    #    biopython: $(python -c "import Bio; print(Bio.__version__)")
    #END_VERSIONS
    '''
}