#!/usr/bin/env nextflow

process DELETE_TEMP_FILES {
    errorStrategy = 'ignore'
    tag "$output_dir"

    input:
    val ready_to_delete
    val directory

    output: 
    val "done", emit: completion_signal

    when: ready_to_delete

    script:
    """
    echo ${directory}
    find "${directory}" -mindepth 1 -maxdepth 1 -type d \\
    ! -name "final_files" \\
    ! -name "failed_assembly" \\
    | xargs -r rm -r \\

    echo "files cleaned!"
    """
}

