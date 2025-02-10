process BWA_MEM2_ALIGN {
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_high'
    container 'ilepeli/revica-strm:0.0.4'

    input:
    tuple val(meta), path(fastq)
    tuple val(meta), val(ref_info), path(ref)

    output:
    tuple val(meta), val(ref_info), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    tuple val(meta), val(ref_info), path(ref),                                      emit: ref
    tuple val(meta), path(fastq),                                                   emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    def input = meta.single_end ? "${fastq}" : "${fastq[0]} ${fastq[1]}"

    // this currently only supports ref paths to a fasta file, not a directory or subdirectory
    """
    bwa-mem2 index $ref

    ## run bwa-mem2 
    bwa-mem2 mem \
        $ref \
        $input \
        -t $task.cpus \
        > ${prefix}.bam

    # remove unmapped reads, sort and index bam files  
    samtools view -b -F 4 -@ ${task.cpus} ${prefix}.bam -o ${prefix}_mapped.bam           
    samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam ${prefix}_mapped.bam  
    samtools index -@ ${task.cpus} ${prefix}.sorted.bam 
    """
}
