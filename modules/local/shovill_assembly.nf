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
    #R1_realpath=$(realpath "!{reads[0]}")
    #R2_realpath=$(realpath "!{reads[0]}")

    shovill -R1 shovill/r1_!{sample}.fastq.gz -R2 shovill/r2_!{sample}.fastq.gz --outdir shovill --force --cpu 4 \
            && mv shovill/contigs.fa !{sample}_assembled.fasta 
            
    rm -rf shovill/r1_!{sample}.fastq.gz shovill/r2_!{sample}.fastq.gz

    '''

}

// process cutadapt {
//     publishDir "${params.final_outdir}/${sample}", mode: 'copy', pattern: "*.fastq"
//     tag "${sample}"
//     cpus = "${params.maxcpus}"
//     memory = "${params.medmems}"
//     errorStrategy 'retry'
//     maxRetries 3
//     // debug true

//     maxForks = "${params.maxcutadapts}"

//     input:
//     tuple val(sample), path(reads)
    
//     output:
//     tuple val(sample), path ("cutadapt/${sample}.trimmed.1.fastq"), path ("cutadapt/${sample}.trimmed.2.fastq"), emit: cutadapt_fastq, optional: true
//     tuple val(sample), path ("cutadapt/${sample}.cutadapt.json"), emit: cutadapt_json, optional: true

//     shell:
//     '''
//     mkdir -p cutadapt
//     run_cutadapt.py -f !{reads[0]} -r !{reads[1]} \
//                     -o cutadapt -s !{sample} -p !{params.primer} \
//                     -e !{params.cutadapt_maxerror} -l !{params.cutadapt_minlength} \
//                     -t !{params.cutadapt_thread} -b !{params.cutadapt_long}

//     '''

// }
