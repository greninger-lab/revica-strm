process CONCAT_INTRASAMPLE_FILES {
    container 'quay.io/biocontainers/mulled-v2-77320db00eefbbf8c599692102c3d387a37ef02a:08144a66f00dc7684fad061f1466033c0176e7ad-0'
    tag "$output_dir"
    label 'process_single'
 
    input: 
    val output_dir
    val ready_to_concat

    when:
    ready_to_concat

    output: 
    path "*_final.fa", emit: merged_fastas


    script: 
    """
    pip install pysam 

    merged_consensus.py $output_dir 

    find $output_dir -type f -name '*_final.fa' -exec mv {} . \\;
    """

}
