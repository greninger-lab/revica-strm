process BAM_TO_FASTQ {

    container 'quay.io/epil02/revica-strm:0.0.4'
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_high'

    input:
    tuple val(meta), val(ref_info), path(bam), path(bai)

    output:
    path "*_SRA.fastq.gz", optional: true, emit: fastq

    script:
    def prefix = task.ext.prefix ?: ''

    """
    if [[ \$(basename "$bam") = "FAILED.sorted.bam" ]]; then
        echo "Skipping bam to fastq conversion; alignment with ${prefix} failed depth/coverage previously"
        exit 0 # shouldn't cause fail if the outputs are optional
    fi

    samtools fastq \\
    $bam \\
    -F 4 \\
    -o "${prefix}_SRA.fastq.gz" \\
    -N \\
    -@ ${task.cpus}
    """
}
