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


    output: 
    path "versions.yml", emit: versions

    script: 
    """
    merged_consensus.py \\
    $output_dir 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
