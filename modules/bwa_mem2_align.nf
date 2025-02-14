process BWA_MEM2_ALIGN {
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_high'
    container 'quay.io/epil02/revica-strm:0.0.4'

    input:
    tuple val(meta), path(fastq)
    tuple val(meta), val(ref_info), path(ref)

    output:
    tuple val(meta), val(ref_info), path("*.sorted.bam"), path("*.sorted.bam.bai"), optional: true, emit: bam
    tuple val(meta), val(ref_info), path(ref),                                      emit: ref
    tuple val(meta), path(fastq),                                                   emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    def input = meta.single_end ? "${fastq}" : "${fastq[0]} ${fastq[1]}"

    def min_coverage = task.ext.min_coverage
    def min_depth = task.ext.min_depth

    // this currently only supports ref paths to a fasta file, not a directory or subdirectory
    """
    bwa-mem2 index $ref

    ## run bwa-mem2 
    bwa-mem2 mem \
        $ref \
        $input \
        -t $task.cpus \
        | samtools view -bS -@ $task.cpus > ${prefix}.bam

    # check if alignment is above depth/coverage thresholds
    pandepth -i ${prefix}.bam -o ${prefix} -t ${task.cpus}
    gunzip ${prefix}.chr.stat.gz 

    prep_pandepth_output.py ${prefix}.chr.stat $ref ${prefix}_covstats.tsv

    coverage=\$(awk 'BEGIN {FS="\t"} NR>1 {print \$5}' "${prefix}_covstats.tsv")
    mean_depth=\$(awk 'BEGIN {FS="\t"} NR>1 {print \$4}' "${prefix}_covstats.tsv")

    # Check if thresholds are met
    if (( \$(echo "\$coverage >= $min_coverage" | bc) == 1 && \$(echo "\$mean_depth >= $min_depth" | bc) == 1 )); then
        echo "alignment thresholds met!"
        samtools view -b -F 4 -@ ${task.cpus} "${prefix}.bam" -o "${prefix}_mapped.bam"
        samtools sort -@ ${task.cpus} -o "${prefix}.sorted.bam" "${prefix}_mapped.bam"
        samtools index -@ ${task.cpus} "${prefix}.sorted.bam"
    else
        echo "alignment failed to meet thresholds! Depth: \$mean_depth, Coverage: \$coverage"
        rm *.bam
        mv ${prefix}_covstats.tsv ${prefix}_failed_assembly.tsv
    fi
    """
}
