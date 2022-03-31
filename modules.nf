process Trimming_SE {
    errorStrategy 'retry'
    maxRetries 1

    input:
        file R1// from input_read_ch
        file ADAPTERS
        val MINLEN
        val SETTING
        val LEADING
        val TRAILING
        val SWINDOW

    output:
        tuple env(base), file("*.trimmed.fastq.gz") // into Trim_out_SE, Trim_out_SE_FQC

    publishDir "${params.outdir}/trimmed_fastqs", mode: 'copy', pattern:'*.trimmed.fastq*'

    script:
    """
    #!/bin/bash
    base=\$(basename ${R1} ".fastq.gz")
    trimmomatic SE -threads ${task.cpus} ${R1} \$base.trimmed.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:${SETTING} LEADING:${LEADING} TRAILING:${TRAILING} SLIDINGWINDOW:${SWINDOW} MINLEN:${MINLEN} 
    """
}

process Aligning_SE {
    errorStrategy 'retry'
    maxRetries 1

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz") // from Trim_out_SE
        file Reference_rv
 	file Reference_hcov
	file Reference_hmpv
	file Reference_hrsv
	file Reference_hpiv

    output:
	tuple val(base), file("${base}_map_all_bbmap_covstats.txt")
	tuple file("rv_ids.txt"), file("hcov_ids.txt"), file("hmpv_ids.txt"), file("hrsv_ids.txt"), file("hpiv_ids.txt")

    publishDir "${params.outdir}/bbmap_covstats_all_ref", mode: 'copy', pattern:'*_map_all_bbmap_covstats.txt'
    
    script:

    """
    #!/bin/bash

    cat ${Reference_rv} ${Reference_hcov} ${Reference_hmpv} ${Reference_hrsv} ${Reference_hpiv} > ${base}_all_ref.fa

    # Isolate accessions of references
    grep ">" ${Reference_rv} | tr -d ">" | cut -d " " -f1 > rv_ids.txt
    grep ">" ${Reference_hcov} | tr -d ">" | cut -d " " -f1 > hcov_ids.txt
    grep ">" ${Reference_hmpv} | tr -d ">" | cut -d " " -f1 > hmpv_ids.txt
    grep ">" ${Reference_hrsv} | tr -d ">" | cut -d " " -f1 > hrsv_ids.txt
    grep ">" ${Reference_hpiv} | tr -d ">" | cut -d " " -f1 > hpiv_ids.txt

    # Map to the multifasta references of rv, hcov, hmpv, hrsv, hpiv
    bbmap.sh \\
        in=${base}.trimmed.fastq.gz \\
        outm=${base}_map_all.sam \\
        ref=${base}_all_ref.fa \\
        threads=${task.cpus} \\
        covstats=${base}_map_all_bbmap_covstats.txt \\
        local=true interleaved=false maxindel=${params.default_maxindel} -Xmx${task.memory.giga}g > ${base}_map_all_stats.txt 2>&1

    # sort bbmap_out.txt based on median_fold (column 10) in descending order 
    head -1 ${base}_map_all_bbmap_covstats.txt > ${base}_map_all_bbmap_covstats.temp.txt
    awk 'NR>1' < ${base}_map_all_bbmap_covstats.txt | sort -nrk10 >> ${base}_map_all_bbmap_covstats.temp.txt
    mv ${base}_map_all_bbmap_covstats.temp.txt ${base}_map_all_bbmap_covstats.txt
    """
}

process Viral_Identification {
    errorStrategy 'retry'
    maxRetries 1

    input:
	tuple val(base), file("${base}_map_all_bbmap_covstats.txt")
	tuple file("rv_ids.txt"), file("hcov_ids.txt"), file("hmpv_ids.txt"), file("hrsv_ids.txt"), file("hpiv_ids.txt")
	
    output:
        file("*vid.txt") optional true
	file("*_failed_assembly.txt") optional true

    publishDir "${params.outdir}/viral_identification", mode: 'copy', pattern:'*vid.txt'
    publishDir "${params.outdir}/failed_assembly", mode: 'copy', pattern:'*_failed_assembly.txt'

    script:
    """
    # find viral reference (-m minimum median_fold stats in bbmap covstats output for a virus to be identified, default 5)
    python3 $workflow.projectDir/bin/viral_identification.py \\
	-bbmap_covstats ${base}_map_all_bbmap_covstats.txt \\
	-b ${base} \\
	-rv_id rv_ids.txt \\
	-hcov_id hcov_ids.txt \\
	-hmpv_id hmpv_ids.txt \\
	-hrsv_id hrsv_ids.txt \\
	-hpiv_id hpiv_ids.txt \\
	-m 5 
    """
}

process Consensus_Generation_Prep_SE {
    errorStrategy 'retry'
    maxRetries 1

    input:
	file "*_vid.txt"
        file Reference_rv
        file Reference_hcov
        file Reference_hmpv
        file Reference_hrsv
        file Reference_hpiv

    output:
	tuple env(base), env(ref_id), env(ref_sp), file("*_*_*.fa") 

    publishDir "${params.outdir}/ref_genome", mode: 'copy', pattern: '*_*_*.fa'

    script:
    """
    base=\$(head -1 *_vid.txt | awk -F '\t' '{print \$1}')
    ref_id=\$(head -1 *_vid.txt | awk -F '\t' '{print \$2}')
    ref_sp=\$(head -1 *_vid.txt | awk -F '\t' '{print \$3}')

    cat ${Reference_rv} ${Reference_hcov} ${Reference_hmpv} ${Reference_hrsv} ${Reference_hpiv} > all_ref.fa
    samtools faidx all_ref.fa \$ref_id > \$base'_'\$ref_id'_'\$ref_sp'.fa'     
    """
}

process Consensus_Generation_SE {
    errorStrategy 'retry'
    maxRetries 1 

    input:
	tuple val(base), val(ref_id), val(ref_sp), file("${base}_${ref_id}_${ref_sp}.fa") 
    output:
	tuple val(base), val(ref_id), val(ref_sp), file("${base}_${ref_id}_${ref_sp}.consensus_final.fa")
	tuple file("*.txt"), file("*.fa"), file("*bam*") 

    publishDir "${params.outdir}/map_ref_stats", mode: 'copy', pattern: '*_map_ref_stats.txt'
    publishDir "${params.outdir}/map_ref_bam_sorted", mode: 'copy', pattern: '*_map_ref.sorted.bam'
    publishDir "${params.outdir}/consensus_final", mode: 'copy', pattern:'*.consensus_final.fa' 
    publishDir "${params.outdir}/consensus_final_bam_sorted", mode: 'copy', pattern:'*_mapf.sorted.bam' 

    script:
    """
    # find real absolute path of outdir
    process_work_dir=\$PWD
    cd $workflow.launchDir
    outdir_realpath=\$(realpath ${params.outdir})
    cd \$process_work_dir

    # get trimmed reads
    cp \${outdir_realpath}/trimmed_fastqs/${base}.trimmed.fastq.gz ${base}.trimmed.fastq.gz

    if [[ ${ref_sp} == rv ]]
    then
    maxindel_num=${params.rv_maxindel}

    elif [[ ${ref_sp} == hcov ]] 
    then
    maxindel_num=${params.hcov_maxindel}

    elif [[ ${ref_sp} == hmpv ]]
    then
    maxindel_num=${params.hmpv_maxindel}

    elif [[ ${ref_sp} == hrsv ]]
    then
    maxindel_num=${params.hrsv_maxindel}

    else
    maxindel_num=${params.default_maxindel}
    fi

    # Map the reads to the reference
    bbmap.sh \\
	in=${base}.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_sp}_map_ref.sam \\
	ref=${base}_${ref_id}_${ref_sp}.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=\$maxindel_num -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_sp}_map_ref_stats.txt 2>&1

    # Convert the output sam file to bam file, sort and index the bam file
    samtools view -S -b -F 4 ${base}_${ref_id}_${ref_sp}_map_ref.sam | samtools sort -@ ${task.cpus} - > ${base}_${ref_id}_${ref_sp}_map_ref.sorted.bam
    rm ${base}_${ref_id}_${ref_sp}_map_ref.sam

    if [[ ${params.deduplicate} == true ]]
    then
    picard MarkDuplicates -I ${base}_${ref_id}_${ref_sp}_map_ref.sorted.bam -O ${base}_${ref_id}_${ref_sp}_map_ref_deduplicated.sorted.bam -M ${base}_${ref_id}_${ref_sp}_picard_output.txt -REMOVE_DUPLICATES true
    mv ${base}_${ref_id}_${ref_sp}_map_ref.sorted.bam ${base}_${ref_id}_${ref_sp}_map_ref_og.sorted.bam
    mv ${base}_${ref_id}_${ref_sp}_map_ref_deduplicated.sorted.bam ${base}_${ref_id}_${ref_sp}_map_ref.sorted.bam
    fi

    # Calling Consensus
    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}_${ref_id}_${ref_sp}.fa \\
        --min-BQ 15 \\
        --output ${base}_${ref_id}_${ref_sp}_1.mpileup \\
        ${base}_${ref_id}_${ref_sp}_map_ref.sorted.bam
    cat ${base}_${ref_id}_${ref_sp}_1.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}_${ref_id}_${ref_sp}.consensus1

    # Get rid of leading and trailing repeated Ns and ns
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}_${ref_id}_${ref_sp}.consensus1.fa > ${base}_${ref_id}_${ref_sp}.consensus1.temp.fa
    cp ${base}_${ref_id}_${ref_sp}.consensus1.temp.fa ${base}_${ref_id}_${ref_sp}.consensus1.fa
    rm ${base}_${ref_id}_${ref_sp}.consensus1.temp.fa

    # Modify first consensus header
    cp ${base}_${ref_id}_${ref_sp}.consensus1.fa ${base}_${ref_id}_${ref_sp}.consensus1.fa.backup
    echo '>${base}_${ref_id}_${ref_sp}.consensus1' > ${base}_${ref_id}_${ref_sp}.consensus1.fa
    tail -n+2 ${base}_${ref_id}_${ref_sp}.consensus1.fa.backup >> ${base}_${ref_id}_${ref_sp}.consensus1.fa
    rm ${base}_${ref_id}_${ref_sp}.consensus1.fa.backup
 
    # Map reads to consensus1 and create bam and sorted bam files
    bbmap.sh \\
	in=${base}.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_sp}_map1.sam \\
	ref=${base}_${ref_id}_${ref_sp}.consensus1.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=\$maxindel_num -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_sp}_map1_stats.txt 2>&1
    samtools view -S -b -F 4 ${base}_${ref_id}_${ref_sp}_map1.sam | samtools sort -@ ${task.cpus} - > ${base}_${ref_id}_${ref_sp}_map1.sorted.bam
    rm ${base}_${ref_id}_${ref_sp}_map1.sam

    # second consensus generation 
    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}_${ref_id}_${ref_sp}.consensus1.fa \\
        --min-BQ 15 \\
        --output ${base}_${ref_id}_${ref_sp}_2.mpileup \\
        ${base}_${ref_id}_${ref_sp}_map1.sorted.bam
    cat ${base}_${ref_id}_${ref_sp}_2.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}_${ref_id}_${ref_sp}.consensus2

    # Get rid of repeated Ns and ns
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}_${ref_id}_${ref_sp}.consensus2.fa > ${base}_${ref_id}_${ref_sp}.consensus2.temp.fa
    cp ${base}_${ref_id}_${ref_sp}.consensus2.temp.fa ${base}_${ref_id}_${ref_sp}.consensus2.fa
    rm ${base}_${ref_id}_${ref_sp}.consensus2.temp.fa

    # Modify second consensus header
    cp ${base}_${ref_id}_${ref_sp}.consensus2.fa ${base}_${ref_id}_${ref_sp}.consensus2.fa.backup
    echo '>${base}_${ref_id}_${ref_sp}.consensus2' > ${base}_${ref_id}_${ref_sp}.consensus2.fa
    tail -n+2 ${base}_${ref_id}_${ref_sp}.consensus2.fa.backup >> ${base}_${ref_id}_${ref_sp}.consensus2.fa
    rm ${base}_${ref_id}_${ref_sp}.consensus2.fa.backup
    
    # Align reads to final consensus and create bam and sorted bam files 
    bbmap.sh \\
	in=${base}.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_sp}_map2.sam \\
	ref=${base}_${ref_id}_${ref_sp}.consensus2.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=\$maxindel_num -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_sp}_map2_stats.txt 2>&1
    samtools view -S -b -F 4 ${base}_${ref_id}_${ref_sp}_map2.sam | samtools sort -@ ${task.cpus} - > ${base}_${ref_id}_${ref_sp}_map2.sorted.bam
    rm ${base}_${ref_id}_${ref_sp}_map2.sam

    # Final Consensus Generation
    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}_${ref_id}_${ref_sp}.consensus2.fa \\
        --min-BQ 15 \\
        --output ${base}_${ref_id}_${ref_sp}_final.mpileup \\
        ${base}_${ref_id}_${ref_sp}_map2.sorted.bam
    cat ${base}_${ref_id}_${ref_sp}_final.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}_${ref_id}_${ref_sp}.consensus_final

    # Get rid of repeated Ns and ns
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}_${ref_id}_${ref_sp}.consensus_final.fa > ${base}_${ref_id}_${ref_sp}.consensus_final.temp.fa
    cp ${base}_${ref_id}_${ref_sp}.consensus_final.temp.fa ${base}_${ref_id}_${ref_sp}.consensus_final.fa
    rm ${base}_${ref_id}_${ref_sp}.consensus_final.temp.fa

    # Modify final consensus header
    cp ${base}_${ref_id}_${ref_sp}.consensus_final.fa ${base}_${ref_id}_${ref_sp}.consensus_final.fa.backup
    echo '>${base}_${ref_id}_${ref_sp}' > ${base}_${ref_id}_${ref_sp}.consensus_final.fa
    tail -n+2 ${base}_${ref_id}_${ref_sp}.consensus_final.fa.backup >> ${base}_${ref_id}_${ref_sp}.consensus_final.fa
    rm ${base}_${ref_id}_${ref_sp}.consensus_final.fa.backup
 
    # Align reads to final consensus and create bam and sorted bam files 
    bbmap.sh \\
	in=${base}.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_sp}_mapf.sam \\
	ref=${base}_${ref_id}_${ref_sp}.consensus_final.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=\$maxindel_num -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_sp}_mapf_stats.txt 2>&1
    samtools view -S -b -F 4 ${base}_${ref_id}_${ref_sp}_mapf.sam | samtools sort -@ ${task.cpus} - > ${base}_${ref_id}_${ref_sp}_mapf.sorted.bam
    rm ${base}_${ref_id}_${ref_sp}_mapf.sam
    """
}

process Serotyping {
    container 'bschiffthaler/ncbi-blast:latest'
    errorStrategy 'retry'
    maxRetries 1

    input:
	tuple val(base), val(ref_id), val(ref_sp), file("${base}_${ref_id}_${ref_sp}.consensus_final.fa")
        file BLASTDB_ALL_1
        file BLASTDB_ALL_2
        file BLASTDB_ALL_3
        file BLASTDB_ALL_4
        file BLASTDB_ALL_5
        file BLASTDB_ALL_6
        file BLASTDB_ALL_7
        file BLASTDB_ALL_8

    output:
        tuple val(base), val(ref_id), val(ref_sp)
        file("*.txt")

    publishDir "${params.outdir}/serotype", mode: 'copy', pattern:'*.txt'

    script:

    """
    #!/bin/bash

    blastn -out ${base}_${ref_id}_${ref_sp}_blast_output.txt -query ${base}_${ref_id}_${ref_sp}.consensus_final.fa -db ${BLASTDB_ALL_1} -outfmt 6 -task blastn -max_target_seqs 1 -evalue 1e-5
    # outfmt6 default values: 'qaccver saccver pident length mismatch gapopen qstart qend sstart send' 

    serotype=\$(awk 'FNR==1{print val,\$2}' ${base}_${ref_id}_${ref_sp}_blast_output.txt | cut -d "_" -f2- | cut -d "/" -f2-)
    echo \$serotype > ${base}_${ref_id}_${ref_sp}_serotype.txt
    """
}

process Summary_Generation {
    errorStrategy 'ignore'

    input:
	tuple val(base), val(ref_id), val(ref_sp)

    output:
	file("${base}_${ref_id}_${ref_sp}_summary.csv")

    script:
    """
    #!/bin/bash
 
    # find real absolute path of outdir
    process_work_dir=\$PWD
    cd $workflow.launchDir
    input_realpath=\$(realpath ${params.reads})
    outdir_realpath=\$(realpath ${params.outdir})
    cd \$process_work_dir
 
    # summary header
    echo Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Species,Mapped_Reads,%_Reads_On_Target,Reference_Length,Percent_Ref_Coverage,Median_Fold,Consensus_Length,Num_Ns,Percent_N,Serotype > ${base}_${ref_id}_${ref_sp}_summary.csv

    # get the number of total reads, trimmed reads, and the percent trimmed
    num_untrimmed=\$((\$(gunzip -c \${input_realpath}/${base}*.fastq.gz | wc -l)/4))
    num_trimmed=\$((\$(gunzip -c \${outdir_realpath}/trimmed_fastqs/${base}*trimmed.fastq.gz | wc -l)/4))
    percent_trimmed=\$((100-\$((100*num_trimmed/num_untrimmed))))

    # get the reference genome (map_all) id, size, coverage
    ref_length=\$(grep ${ref_id} \${outdir_realpath}/bbmap_covstats_all_ref/${base}_map_all_bbmap_covstats.txt | awk -F '\t' '{print \$3}')
    ref_coverage=\$(grep ${ref_id} \${outdir_realpath}/bbmap_covstats_all_ref/${base}_map_all_bbmap_covstats.txt | awk -F '\t' '{print \$5}')
    median_fold=\$(grep ${ref_id} \${outdir_realpath}/bbmap_covstats_all_ref/${base}_map_all_bbmap_covstats.txt | awk -F '\t' '{print \$10}')

    # get the number and the percentage of trimmed reads mapped to the reference genome (map_ref)
    mapped_reads=\$(cat \${outdir_realpath}/map_ref_stats/${base}_${ref_id}_${ref_sp}_map_ref_stats.txt | grep "mapped:" | cut -d\$'\\t' -f3)
    mapped_reads=\$(echo \$mapped_reads | awk '{print \$1+\$2}')
    percent_mapped_reads=\$(echo "\$mapped_reads/\$num_trimmed*100" | bc -l | awk 'FNR==1{print val,\$1}')

    # get the consensus final length
    consensus_length=\$(awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length(\$0)}END{print l}' \${outdir_realpath}/consensus_final/${base}_${ref_id}_${ref_sp}.consensus_final.fa | awk 'FNR==2{print val,\$1}')

    # get the number of Ns in the consensus final 
    num_ns=\$(grep -v "^>" \${outdir_realpath}/consensus_final/${base}_${ref_id}_${ref_sp}.consensus_final.fa | tr -cd N | wc -c | awk 'FNR==1{print val,\$1}')

    # get the percentage of Ns in the consensus final
    percent_n=\$(echo "\$num_ns/\$consensus_length*100" | bc -l | awk 'FNR==1{print val,\$1}')

    # get the serotype
    serotype=\$(awk 'FNR==1{print \$1}' \${outdir_realpath}/serotype/${base}_${ref_id}_${ref_sp}_serotype.txt)

    printf "${base},\$num_untrimmed,\$num_trimmed,\$percent_trimmed,${ref_id},${ref_sp},\$mapped_reads,\$percent_mapped_reads,\$ref_length,\$ref_coverage,\$median_fold,\$consensus_length,\$num_ns,\$percent_n,\$serotype" >> ${base}_${ref_id}_${ref_sp}_summary.csv
    """
}

process Final_Processing {
    errorStrategy 'ignore'

    input:
	file("*_summary.csv")

    output: 
	file("summary.csv")

    publishDir "${params.outdir}/run_summary", mode: 'copy', pattern:'summary.csv'

    script:
    """
    #!/bin/bash

    echo Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Species,Mapped_Reads,%_Reads_On_Target,Reference_Length,Percent_Ref_Coverage,Median_Fold,Consensus_Length,Num_Ns,Percent_N,Serotype > summary.csv

    awk '(NR == 2) || (FNR > 1)' *_summary.csv >> summary.csv 

    """
}

process Trimming_PE { 
    errorStrategy 'retry'
    maxRetries 1

    input:
        tuple val(base), file(R1), file(R2) // from input_read_ch
        file ADAPTERS
        val MINLEN
        val SETTING
        val LEADING
        val TRAILING
        val SWINDOW
    output: 
        tuple val(base), file("${base}.R1.paired.trimmed.fastq.gz"), file("${base}.R2.paired.trimmed.fastq.gz")

    publishDir "${params.outdir}/trimmed_fastqs", mode: 'copy',pattern:'*.paired.trimmed.fastq.gz'

    script:
    """
    #!/bin/bash
    trimmomatic PE -threads ${task.cpus} ${R1} ${R2} ${base}.R1.paired.trimmed.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}.R2.paired.trimmed.fastq.gz ${base}.R2.unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:${SETTING} LEADING:${LEADING} TRAILING:${TRAILING} SLIDINGWINDOW:${SWINDOW} MINLEN:${MINLEN}
    """
}

process Aligning_PE {
    errorStrategy 'retry'
    maxRetries 1

    input: 
        tuple val(base), file("${base}.R1.paired.trimmed.fastq.gz"), file("${base}.R2.paired.trimmed.fastq.gz")
        file Reference_rv
        file Reference_hcov
	file Reference_hmpv
	file Reference_hrsv
	file Reference_hpiv

    output:
	tuple val(base), file("${base}_map_all_bbmap_covstats.txt")
	tuple file("rv_ids.txt"), file("hcov_ids.txt"), file("hmpv_ids.txt"), file("hrsv_ids.txt"), file("hpiv_ids.txt")

    publishDir "${params.outdir}/bbmap_covstats_all_ref", mode: 'copy', pattern:'*_map_all_bbmap_covstats.txt'
    
    script:

    """
    #!/bin/bash

    cat ${Reference_rv} ${Reference_hcov} ${Reference_hmpv} ${Reference_hrsv} ${Reference_hpiv} > ${base}_all_ref.fa

    # Isolate accessions of references
    grep ">" ${Reference_rv} | tr -d ">" | cut -d " " -f1 > rv_ids.txt
    grep ">" ${Reference_hcov} | tr -d ">" | cut -d " " -f1 > hcov_ids.txt
    grep ">" ${Reference_hmpv} | tr -d ">" | cut -d " " -f1 > hmpv_ids.txt
    grep ">" ${Reference_hrsv} | tr -d ">" | cut -d " " -f1 > hrsv_ids.txt
    grep ">" ${Reference_hpiv} | tr -d ">" | cut -d " " -f1 > hpiv_ids.txt

    # Map to the multifasta references of rv, hcov, hmpv, hrsv, hpiv
    bbmap.sh \\
        in=${base}.R1.paired.trimmed.fastq.gz \\
	in2=${base}.R2.paired.trimmed.fastq.gz \\
        outm=${base}_map_all.sam \\
        ref=${base}_all_ref.fa \\
        threads=${task.cpus} \\
        covstats=${base}_map_all_bbmap_covstats.txt \\
        local=true interleaved=false maxindel=80 -Xmx${task.memory.giga}g > ${base}_map_all_stats.txt 2>&1

    # sort bbmap_out.txt based on median_fold (column 10) in descending order 
    head -1 ${base}_map_all_bbmap_covstats.txt > ${base}_map_all_bbmap_covstats.temp.txt
    awk 'NR>1' < ${base}_map_all_bbmap_covstats.txt | sort -nrk10 >> ${base}_map_all_bbmap_covstats.temp.txt
    mv ${base}_map_all_bbmap_covstats.temp.txt ${base}_map_all_bbmap_covstats.txt
    """
}

process Consensus_Generation_Prep_PE {
    errorStrategy 'retry'
    maxRetries 1

    input:
	file "*_vid.txt"
        file Reference_rv
        file Reference_hcov
        file Reference_hmpv
        file Reference_hrsv
	file Reference_hpiv

    output:
	tuple env(base), env(ref_id), env(ref_sp), file("*_*_*.fa") 

    publishDir "${params.outdir}/ref_genome", mode: 'copy', pattern: '*_*_*.fa'

    script:
    """
    base=\$(head -1 *_vid.txt | awk '{print \$1}')
    ref_id=\$(head -1 *_vid.txt | awk '{print \$2}')
    ref_sp=\$(head -1 *_vid.txt | awk '{print \$3}')

    cat ${Reference_rv} ${Reference_hcov} ${Reference_hmpv} ${Reference_hrsv} ${Reference_hpiv} > all_ref.fa
    samtools faidx all_ref.fa \$ref_id > \$base'_'\$ref_id'_'\$ref_sp'.fa'     
    """
}

process Consensus_Generation_PE {
    errorStrategy 'retry'
    maxRetries 1 

    input:
	tuple val(base), val(ref_id), val(ref_sp), file("${base}_${ref_id}_${ref_sp}.fa") 
    output:
	tuple val(base), val(ref_id), val(ref_sp), file("${base}_${ref_id}_${ref_sp}.consensus_final.fa")
	tuple file("*.txt"), file("*.mpileup"), file("*.fa"), file("*bam*") 

    publishDir "${params.outdir}/map_ref_stats", mode: 'copy', pattern: '*_map_ref_stats.txt'
    publishDir "${params.outdir}/map_ref_bam_sorted", mode: 'copy', pattern: '*_map_ref.sorted.bam'
    publishDir "${params.outdir}/consensus_final", mode: 'copy', pattern:'*.consensus_final.fa' 
    publishDir "${params.outdir}/consensus_final_bam_sorted", mode: 'copy', pattern:'*_mapf.sorted.bam' 

    script:
    """

    # find real absolute path of outdir
    process_work_dir=\$PWD
    cd $workflow.launchDir
    outdir_realpath=\$(realpath ${params.outdir})
    cd \$process_work_dir
    
    # get paired trimmed reads
    cp \${outdir_realpath}/trimmed_fastqs/${base}.R1.paired.trimmed.fastq.gz ${base}.R1.paired.trimmed.fastq.gz
    cp \${outdir_realpath}/trimmed_fastqs/${base}.R2.paired.trimmed.fastq.gz ${base}.R2.paired.trimmed.fastq.gz

    if [[ ${ref_sp} == rv ]]
    then
    maxindel_num=${params.rv_maxindel}

    elif [[ ${ref_sp} == hcov ]] 
    then
    maxindel_num=${params.hcov_maxindel}

    elif [[ ${ref_sp} == hmpv ]]
    then
    maxindel_num=${params.hmpv_maxindel}

    elif [[ ${ref_sp} == hrsv ]]
    then
    maxindel_num=${params.hrsv_maxindel}

    else
    maxindel_num=${params.default_maxindel}
    fi

    # Map the reads to the reference
    bbmap.sh \\
	in=${base}.R1.paired.trimmed.fastq.gz \\
	in2=${base}.R2.paired.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_sp}_map_ref.sam \\
	ref=${base}_${ref_id}_${ref_sp}.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=\$maxindel_num -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_sp}_map_ref_stats.txt 2>&1

    # Convert the output sam file to bam file, sort and index the bam file
    samtools view -S -b -F 4 ${base}_${ref_id}_${ref_sp}_map_ref.sam | samtools sort -@ ${task.cpus} - > ${base}_${ref_id}_${ref_sp}_map_ref.sorted.bam
    rm ${base}_${ref_id}_${ref_sp}_map_ref.sam

    if [[ ${params.deduplicate} == true ]]
    then
    picard MarkDuplicates -I ${base}_${ref_id}_${ref_sp}_map_ref.sorted.bam -O ${base}_${ref_id}_${ref_sp}_map_ref_deduplicated.sorted.bam -M ${base}_${ref_id}_${ref_sp}_picard_output.txt -REMOVE_DUPLICATES true
    mv ${base}_${ref_id}_${ref_sp}_map_ref.sorted.bam ${base}_${ref_id}_${ref_sp}_map_ref_og.sorted.bam
    mv ${base}_${ref_id}_${ref_sp}_map_ref_deduplicated.sorted.bam ${base}_${ref_id}_${ref_sp}_map_ref.sorted.bam
    fi

    # Calling Consensus
    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}_${ref_id}_${ref_sp}.fa \\
        --min-BQ 15 \\
        --output ${base}_${ref_id}_${ref_sp}_1.mpileup \\
        ${base}_${ref_id}_${ref_sp}_map_ref.sorted.bam
    cat ${base}_${ref_id}_${ref_sp}_1.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}_${ref_id}_${ref_sp}.consensus1

    # Get rid of leading and trailing repeated Ns and ns
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}_${ref_id}_${ref_sp}.consensus1.fa > ${base}_${ref_id}_${ref_sp}.consensus1.temp.fa
    mv ${base}_${ref_id}_${ref_sp}.consensus1.temp.fa ${base}_${ref_id}_${ref_sp}.consensus1.fa
 
    # Modify first consensus header
    cp ${base}_${ref_id}_${ref_sp}.consensus1.fa ${base}_${ref_id}_${ref_sp}.consensus1.fa.backup
    echo '>${base}_${ref_id}_${ref_sp}.consensus1' > ${base}_${ref_id}_${ref_sp}.consensus1.fa
    tail -n+2 ${base}_${ref_id}_${ref_sp}.consensus1.fa.backup >> ${base}_${ref_id}_${ref_sp}.consensus1.fa
    rm ${base}_${ref_id}_${ref_sp}.consensus1.fa.backup

    # Map reads to consensus1 and create bam and sorted bam files
    bbmap.sh \\
	in=${base}.R1.paired.trimmed.fastq.gz \\
	in2=${base}.R2.paired.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_sp}_map1.sam \\
	ref=${base}_${ref_id}_${ref_sp}.consensus1.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=\$maxindel_num -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_sp}_map1_stats.txt 2>&1
    samtools view -S -b -F 4 ${base}_${ref_id}_${ref_sp}_map1.sam | samtools sort -@ ${task.cpus} - > ${base}_${ref_id}_${ref_sp}_map1.sorted.bam
    rm ${base}_${ref_id}_${ref_sp}_map1.sam

    # second consensus generation 
    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}_${ref_id}_${ref_sp}.consensus1.fa \\
        --min-BQ 15 \\
        --output ${base}_${ref_id}_${ref_sp}_2.mpileup \\
        ${base}_${ref_id}_${ref_sp}_map1.sorted.bam
    cat ${base}_${ref_id}_${ref_sp}_2.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}_${ref_id}_${ref_sp}.consensus2

    # Get rid of repeated Ns and ns
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}_${ref_id}_${ref_sp}.consensus2.fa > ${base}_${ref_id}_${ref_sp}.consensus2.temp.fa
    mv ${base}_${ref_id}_${ref_sp}.consensus2.temp.fa ${base}_${ref_id}_${ref_sp}.consensus2.fa
    
    # Modify second consensus header
    cp ${base}_${ref_id}_${ref_sp}.consensus2.fa ${base}_${ref_id}_${ref_sp}.consensus2.fa.backup
    echo '>${base}_${ref_id}_${ref_sp}.consensus2' > ${base}_${ref_id}_${ref_sp}.consensus2.fa
    tail -n+2 ${base}_${ref_id}_${ref_sp}.consensus2.fa.backup >> ${base}_${ref_id}_${ref_sp}.consensus2.fa
    rm ${base}_${ref_id}_${ref_sp}.consensus2.fa.backup

    # Align reads to final consensus and create bam and sorted bam files 
    bbmap.sh \\
	in=${base}.R1.paired.trimmed.fastq.gz \\
	in2=${base}.R2.paired.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_sp}_map2.sam \\
	ref=${base}_${ref_id}_${ref_sp}.consensus2.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=\$maxindel_num -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_sp}_map2_stats.txt 2>&1
    samtools view -S -b -F 4 ${base}_${ref_id}_${ref_sp}_map2.sam | samtools sort -@ {$task.cpus} - > ${base}_${ref_id}_${ref_sp}_map2.sorted.bam
    rm ${base}_${ref_id}_${ref_sp}_map2.sam

    # Final Consensus Generation
    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}_${ref_id}_${ref_sp}.consensus2.fa \\
        --min-BQ 15 \\
        --output ${base}_${ref_id}_${ref_sp}_final.mpileup \\
        ${base}_${ref_id}_${ref_sp}_map2.sorted.bam
    cat ${base}_${ref_id}_${ref_sp}_final.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}_${ref_id}_${ref_sp}.consensus_final

    # Get rid of repeated Ns and ns
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}_${ref_id}_${ref_sp}.consensus_final.fa > ${base}_${ref_id}_${ref_sp}.consensus_final.temp.fa
    mv ${base}_${ref_id}_${ref_sp}.consensus_final.temp.fa ${base}_${ref_id}_${ref_sp}.consensus_final.fa
 
    # Modify final consensus header
    cp ${base}_${ref_id}_${ref_sp}.consensus_final.fa ${base}_${ref_id}_${ref_sp}.consensus_final.fa.backup
    echo '>${base}_${ref_id}_${ref_sp}' > ${base}_${ref_id}_${ref_sp}.consensus_final.fa
    tail -n+2 ${base}_${ref_id}_${ref_sp}.consensus_final.fa.backup >> ${base}_${ref_id}_${ref_sp}.consensus_final.fa
    rm ${base}_${ref_id}_${ref_sp}.consensus_final.fa.backup

    # Align reads to final consensus and create bam and sorted bam files 
    bbmap.sh \\
	in=${base}.R1.paired.trimmed.fastq.gz \\
	in2=${base}.R2.paired.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_sp}_mapf.sam \\
	ref=${base}_${ref_id}_${ref_sp}.consensus_final.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=\$maxindel_num -Xmx${task.memory.giga}g > ${base}_mapf_stats.txt 2>&1
    samtools view -S -b -F 4 ${base}_${ref_id}_${ref_sp}_mapf.sam | samtools sort -@ ${task.cpus} - > ${base}_${ref_id}_${ref_sp}_mapf.sorted.bam
    rm ${base}_${ref_id}_${ref_sp}_mapf.sam
    """
}
