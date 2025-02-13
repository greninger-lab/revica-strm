process BWA_MEM2_ALIGN_DB {
    tag "$meta.id"
    label 'process_high'
    container 'ilepeli/revica-strm:0.0.4'
    maxForks 2 // Still testing this out


    input:
    tuple val(meta), path(fastq)
    tuple path(db), path(db_indexed)

    output:
    tuple val(meta), path("*covstats.tsv"), emit: covstats

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "$fastq" : "${fastq[0]} ${fastq[1]}"

    """
    ## run bwa-mem2 
    bwa-mem2 mem \
        $db \
        $input \
        -t $task.cpus \
        | samtools view -bS -@ $task.cpus > ${prefix}.bam


    ## run pandepth to get coverage/depth reporting
    pandepth -i ${prefix}.bam -o ${prefix} -t ${task.cpus}
    gunzip ${prefix}.chr.stat.gz

    ## replace abbreviated ref names in pandepth with originals from db
    ## we do this because BWA_MEM only records the alignment ref before the first space, which is
    ## usually just the acc number. We need the rest of the fasta header for downstream analyses

    prep_pandepth_output.py ${prefix}.chr.stat $db ${prefix}_covstats.tsv
    """
}
