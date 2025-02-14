process BWA_MEM2_ALIGN_DB {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/epil02/revica-strm:0.0.4'


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
    // this flag combines 0x4 and 0x800, which means:
    // 1. read is unmapped
    // 2. read is supplementary alignment (chimeric, not representative alignment)
    // we can't let these reads survive.
    FLAG=2052

    ## run bwa-mem2 
    bwa-mem2 mem \
        $db \
        $input \
        -t $task.cpus \
        | samtools view -bS -F \$FLAG -@ 2 > ${prefix}.bam


    ## run pandepth to get coverage/depth reporting
    pandepth -i ${prefix}.bam -o ${prefix} -t ${task.cpus}
    gunzip ${prefix}.chr.stat.gz

    ## replace abbreviated ref names in pandepth with originals from db
    ## we do this because BWA_MEM only records the alignment ref before the first space, which is
    ## usually just the acc number. We need the rest of the fasta header for downstream analyses

    prep_pandepth_output.py ${prefix}.chr.stat $db ${prefix}_covstats.tsv
    """
}
