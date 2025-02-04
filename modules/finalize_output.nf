process FINALIZE_OUTPUT {
    // container 'quay.io/biocontainers/mulled-v2-77320db00eefbbf8c599692102c3d387a37ef02a:08144a66f00dc7684fad061f1466033c0176e7ad-0'
    container 'greningerlab/revica:ubuntu-20.04'
    tag "$output_dir"
    label 'process_medium'
 
    input: 
    val output_dir
    val samplesheet
    val ready_to_concat
    val concat_flu

    when:
    ready_to_concat

    output: 
    path "*.fa",            emit: fastas,   optional: true
    path "*.bam",           emit: bams,     optional: true
    path "*.bai",           emit: bais,     optional: true
    path "*SRA.fastq.gz",   emit: fastqs
    val "done",             emit: done

    script: 
    """
    if [ "${concat_flu}" = "true" ]; then 
        finalize_output.py $output_dir $samplesheet $task.cpus --merge
    else
        finalize_output.py $output_dir $samplesheet $task.cpus
    fi
    """

}
