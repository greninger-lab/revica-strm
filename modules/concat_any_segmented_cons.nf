process CONCAT_ANY_SEGMENTED_CONS {

    debug true
 
    tag "$output_dir"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.8.3' 

    input: 
    val output_dir
    val ready_to_concat

    when:
    ready_to_concat

    script: 
    """
    merged_consensus.py \\
    $output_dir 
    """

}
