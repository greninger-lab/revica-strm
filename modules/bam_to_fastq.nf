process BAM_TO_FASTQ {

    container 'quay.io/epil02/revica-strm:0.0.4'
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_high'

    input:
    tuple val(meta), val(ref_info), path(bam), path(bai)
    val include_unpaired

    output:
    path "*_SRA.fastq.gz", optional: true, emit: fastq

    script:
    def prefix = task.ext.prefix ?: ''

    // if single end or include unpaired, output everything to a single fastq file
    // otherwise, if paired end and not including unpaired, output to -1 and -2.
    def output = !meta.single_end && !include_unpaired ? "-1 ${prefix}_SRA_1.fastq.gz -2 ${prefix}_SRA_1.fastq.gz -s ${prefix}_unpaired.fastq.gz" : " -0 ${prefix}_SRA.fastq.gz -s ${prefix}_SRA.fastq.gz"


    """
    if [[ \$(basename "$bam") = "FAILED.sorted.bam" ]]; then
        echo "Skipping bam to fastq conversion; alignment with ${prefix} failed depth/coverage previously"
        exit 0 # shouldn't cause fail if the outputs are optional
    fi

    samtools fastq \\
    $bam \\
    -F 4 \\
    $output \\
    -N \\
    -@ ${task.cpus}
    """
}
