process BAM_TO_FASTQ {

    container 'greningerlab/revica:ubuntu-20.04'
    tag "$output_dir"
    label 'process_single'

    input:
    path bam

    output:
    path "*_SRA_.fastq.gz", emit: fastq

    script:

    """

    samtools fastq \\
    $bam \\
    -F 4 \\
    -o "${bam.baseName}_SRA_.fastq.gz" \\
    -N \\
    -@ ${task.cpus}
    """
}
