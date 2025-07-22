process BAM_TO_FASTQ {

    container 'quay.io/epil02/revica-strm:0.0.6'
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_high'

    input:
    tuple val(meta), val(ref_info), path(bam), path(bai)
    val include_unpaired

    output:
    path "*_SRA.fastq.gz", optional: true, emit: fastq

    script:
    def prefix = task.ext.prefix ?: ''

    def output = ""
    def prep = ""
    def inbam = ""

    // single end, including everything
    if (meta.single_end) {
        output = "-0 ${prefix}_SRA.fastq.gz"
        prep = ""
        inbam = "$bam"
    }

    // paired end, including everything, change flag
    if (!include_unpaired & !meta.single_end) {
        output = "-1 ${prefix}_1_SRA.fastq.gz -2 ${prefix}_2_SRA.fastq.gz"
        prep = "samtools sort -o ${prefix}_name_sorted.bam -N -@ ${task.cpus} $bam" 
        inbam = "${prefix}_name_sorted.bam"
    }

    // paired end, include pairs only
    if (include_unpaired & !meta.single_end) {
        output = "-1 ${prefix}_1_SRA.fastq.gz -2 ${prefix}_2_SRA.fastq.gz -s ${prefix}_unpaired_SRA.fastq.gz"
        prep = "samtools sort -o ${prefix}_name_sorted.bam -N -@ ${task.cpus} $bam" 
        inbam = "${prefix}_name_sorted.bam"
    }

    """

    echo $output
    echo $include_unpaired

    if [[ \$(basename "$bam") = "FAILED.sorted.bam" ]]; then
        echo "Skipping bam to fastq conversion; alignment with ${prefix} failed depth/coverage previously"
        exit 0 # shouldn't cause fail if the outputs are optional
    fi

    $prep

    samtools fastq \\
    $inbam \\
    -F 4 \\
    -N \\
    -@ ${task.cpus} \\
    $output

    rm -f ${prefix}_name_sorted.bam

    """
}
