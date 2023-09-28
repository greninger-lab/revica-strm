process BBMAP_ALIGN_DB {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/mulled-v2-008daec56b7aaf3f162d7866758142b9f889d690:e8a286b2e789c091bac0a57302cdc78aa0112353-0'

    input:
    tuple val(meta), path(fastq)
    path db

    output:
    tuple val(meta), path("*_covstats.tsv"),    emit: covstats
    tuple val(meta), path("*.log"),             emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "in=${fastq}" : "in=${fastq[0]} in2=${fastq[1]}"

    // Set the ref_seq variable to reflect the three possible types of reference input: 1) directory
    // named 'ref', 2) directory named something else (containg a 'ref' subdir) or 3) a sequence
    // file in fasta format
    if ( db.isDirectory() ) {
        if ( db ==~ /(.\/)?db\/?/ ) {
            ref_seq = ''
        } else {
            ref_seq = "path=${db}"
        }
    } else {
        ref_seq = "ref=${db}"
    }

    """
    bbmap.sh \\
        ${ref_seq} \\
        $input \\
        out=${prefix}.bam \\
        $args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g \\
        covstats=${prefix}_covstats.tsv \\
        &> ${prefix}.bbmap.log

    # sort bbmap_out.txt based on median_fold (column 10) in descending order 
    head -1 ${prefix}_covstats.tsv > ${prefix}_covstats.temp.txt 
    awk 'NR>1' < ${prefix}_covstats.tsv | sort -t \$'\t' -nrk10 >> ${prefix}_covstats.temp.txt
    mv ${prefix}_covstats.temp.txt ${prefix}_covstats.tsv
    """
}
