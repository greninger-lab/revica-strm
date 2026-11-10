process DEDUP_CONSENSUS {
    tag "${meta.id}_${dupchk_tag}"
    label 'process_single'
    container 'quay.io/biocontainers/pywfa:0.6.0--py312hfabe715_0'

    input:
    tuple val(meta), val(dupchk_tag), val(ref_infos), path(covstats), path(fastas)

    output:
    path("*_dedup_report.tsv"), emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def tags  = ref_infos.collect { it.tag }.join(" ")
    def threshold = params.dedup_threshold
    """
    dedup_consensus.py \\
        --fastas ${fastas} \\
        --covstats ${covstats} \\
        --ref-tags ${tags} \\
        --basename ${meta.id} \\
        --dupchk-tag ${dupchk_tag} \\
        --threshold ${threshold}
    """
}
