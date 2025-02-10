process BWA_MEM2_INDEX {
    tag "${fasta}"
    label 'process_single'
    container 'ilepeli/revica-strm:0.0.4'

    input: 
    path fasta

    output:
    tuple path(fasta), path("${fasta}*"), emit: indexed_fasta

    """
    bwa-mem2 index $fasta
    """
}

    
