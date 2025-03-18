process shovill_assembly {
    publishDir "${params.final_outdir}/${sample}", mode: 'copy', pattern: "*.fasta"
    tag "${sample}"
    // debug true

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path ("${sample}_assembled.fasta"), emit: fasta, optional:true


    shell:
    '''
    # comment here
    mkdir -p shovill
    cp !{reads[0]} shovill/r1_!{sample}.fastq.gz
    cp !{reads[1]} shovill/r2_!{sample}.fastq.gz

    shovill -R1 shovill/r1_!{sample}.fastq.gz -R2 shovill/r2_!{sample}.fastq.gz --outdir shovill --force --cpu 4 \
            && mv shovill/contigs.fa !{sample}_assembled.fasta 
            
    rm -rf shovill/r1_!{sample}.fastq.gz shovill/r2_!{sample}.fastq.gz

    '''

}
