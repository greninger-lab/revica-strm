process Trimming_SE {
    // errorStrategy 'ignore'
    errorStrategy 'retry'
    maxRetries 1


    input:
        file R1// from input_read_ch
        file ADAPTERS

    output:
        tuple env(base), file("*.trimmed.fastq.gz") 
	file("*_trim_stats.txt")

    publishDir "${params.outdir}/trimmed_fastq", mode: 'copy', pattern:'*.trimmed.fastq.gz'
    publishDir "${params.outdir}/trim_stats", mode: 'copy', pattern:'*_trim_stats.txt'

    script:
    """
    #!/bin/bash

    if [[ ${R1} == *.fastq ]]
    then
    base=\$(basename ${R1} ".fastq")
    	if [[ ${params.sample} == false ]]
    	then 
            gzip ${R1} 
    	    java -Xmx${task.memory.giga}g -jar \$TRIMMOMATIC SE -threads ${task.cpus} ${R1}.gz \${base}.trimmed.fastq.gz ILLUMINACLIP:${ADAPTERS}:${params.SETTING} \
    	    LEADING:${params.LEADING} TRAILING:${params.TRAILING} SLIDINGWINDOW:${params.SWINDOW} MINLEN:${params.MINLEN} &> \${base}_trim_stats.txt
	else
	    seqtk sample ${R1} ${params.sample} | gzip > \${base}_sampled.fastq.gz
    	    java -Xmx${task.memory.giga}g -jar \$TRIMMOMATIC SE -threads ${task.cpus} \${base}_sampled.fastq.gz \${base}.trimmed.fastq.gz ILLUMINACLIP:${ADAPTERS}:${params.SETTING} \
    	    LEADING:${params.LEADING} TRAILING:${params.TRAILING} SLIDINGWINDOW:${params.SWINDOW} MINLEN:${params.MINLEN} &> \${base}_trim_stats.txt
	fi
    elif [[ ${R1} == *.fastq.gz ]] 
    then
    base=\$(basename ${R1} ".fastq.gz")
	if [[ ${params.sample} == false ]]
	then
    	    java -Xmx${task.memory.giga}g -jar \$TRIMMOMATIC SE -threads ${task.cpus} ${R1} \${base}.trimmed.fastq.gz ILLUMINACLIP:${ADAPTERS}:${params.SETTING} \
    	    LEADING:${params.LEADING} TRAILING:${params.TRAILING} SLIDINGWINDOW:${params.SWINDOW} MINLEN:${params.MINLEN} &> \${base}_trim_stats.txt
	else
	    seqtk sample ${R1} ${params.sample} | gzip > \${base}_sampled.fastq.gz
    	    java -Xmx${task.memory.giga}g -jar \$TRIMMOMATIC SE -threads ${task.cpus} \${base}_sampled.fastq.gz \${base}.trimmed.fastq.gz ILLUMINACLIP:${ADAPTERS}:${params.SETTING} \
    	    LEADING:${params.LEADING} TRAILING:${params.TRAILING} SLIDINGWINDOW:${params.SWINDOW} MINLEN:${params.MINLEN} &> \${base}_trim_stats.txt
	fi
    fi
    """
}

process Aligning_SE {
    errorStrategy 'retry'
    maxRetries 1

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz")
	file ref

    output:
	tuple val(base), file("${base}_map_all_bbmap_covstats.tsv"), file("${base}_map_all_stats.txt")

    publishDir "${params.outdir}/bbmap_covstats_all_ref", mode: 'copy', pattern:'*_map_all_bbmap_covstats.tsv'
    
    script:

    """
    #!/bin/bash

    # Map to the multifasta reference
    bbmap.sh \\
        in=${base}.trimmed.fastq.gz \\
        outm=${base}_map_all.sam \\
        ref=${ref} \\
        threads=${task.cpus} \\
        covstats=${base}_map_all_bbmap_covstats.tsv \\
        local=true interleaved=false maxindel=80 -Xmx${task.memory.giga}g > ${base}_map_all_stats.txt 2>&1

    # sort bbmap_out.txt based on median_fold (column 10) in descending order 
    head -1 ${base}_map_all_bbmap_covstats.tsv > ${base}_map_all_bbmap_covstats.temp.txt
    awk 'NR>1' < ${base}_map_all_bbmap_covstats.tsv | sort -t \$'\t' -nrk10 >> ${base}_map_all_bbmap_covstats.temp.txt
    mv ${base}_map_all_bbmap_covstats.temp.txt ${base}_map_all_bbmap_covstats.tsv

    """
}

process Viral_Identification {
    errorStrategy 'retry'
    maxRetries 1

    input:
	tuple val(base), file("${base}_map_all_bbmap_covstats.tsv"), file("${base}_map_all_stats.txt")
	file ref
	
    output:
        file("*vid.txt") optional true
	file("*_failed_assembly.tsv") optional true

    publishDir "${params.outdir}/viral_identification", mode: 'copy', pattern:'*vid.txt'
    publishDir "${params.outdir}/failed_assembly", mode: 'copy', pattern:'*_failed_assembly.tsv'

    script:
    """
    #!/bin/bash

    # get total input reads for mapping
    reads_used=\$(grep "Reads Used:" ${base}_map_all_stats.txt | cut -f2)

    # get mapped reads
    mapped_reads=\$(grep "Mapped reads:" ${base}_map_all_stats.txt | cut -f2)
    
    # identify initial reference
    python3 $workflow.projectDir/bin/select_reference.py -bbmap_covstats ${base}_map_all_bbmap_covstats.tsv -b ${base} -reads_num \${reads_used} -mapped_reads \${mapped_reads} -m ${params.m} -p ${params.p}

    """
}

process Consensus_Generation_Prep_SE {
    errorStrategy 'retry'
    maxRetries 1

    input:
	file "*_vid.txt"
	file ref

    output:
	tuple env(base), env(ref_id), env(ref_tag), file("*.fa") 

    publishDir "${params.outdir}/ref_genome", mode: 'copy', pattern: '*.fa'

    script:
    """
    #!/bin/bash

    base=\$(head -1 *_vid.txt | cut -f1)
    ref_id=\$(head -1 *_vid.txt | cut -f2)
    ref_tag=\$(head -1 *_vid.txt | cut -f3)
    samtools faidx ${ref} \$ref_id > \${base}'_'\${ref_id}'_'\${ref_tag}'.fa'     

    """
}

process Consensus_Generation_SE {
    errorStrategy 'retry'
    maxRetries 1 

    input:
	tuple val(base), val(ref_id), val(ref_tag), file("${base}_${ref_id}_${ref_tag}.fa") 

    output:
        tuple val(base), val(ref_id), val(ref_tag), file("${base}_${ref_id}_${ref_tag}.fa"), file("${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam")
	tuple val(base), val(ref_id), val(ref_tag), file("${base}_${ref_id}_${ref_tag}.consensus_final.fa")
	tuple file("*.txt"), file("*.fa"), file("*bam*"), file("*.bai") 

    publishDir "${params.outdir}/map_ref_bam_sorted", mode: 'copy', pattern: '*_map_ref.sorted.bam'
    publishDir "${params.outdir}/map_ref_bam_sorted", mode: 'copy', pattern: '*_map_ref.sorted.bam.bai'
    publishDir "${params.outdir}/consensus_final", mode: 'copy', pattern:'*.consensus_final.fa' 
    publishDir "${params.outdir}/consensus_final_bam_sorted", mode: 'copy', pattern:'*_mapf.sorted.bam' 
    publishDir "${params.outdir}/consensus_final_bam_sorted", mode: 'copy', pattern:'*_mapf.sorted.bam.bai' 
    publishDir "${params.outdir}/map_consensus_final_stats", mode: 'copy', pattern: '*_mapf_stats.txt'

    script:
    """
    #!/bin/bash

    # find real absolute path of outdir
    process_work_dir=\$PWD
    cd ${workflow.launchDir}
    outdir_realpath=\$(realpath ${params.outdir})
    cd \$process_work_dir

    # Map the reads to the reference
    bbmap.sh \\
	in=\${outdir_realpath}/trimmed_fastq/${base}.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_tag}_map_ref.sam \\
	ref=${base}_${ref_id}_${ref_tag}.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=80 ambiguous=random -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_tag}_map_ref_stats.txt 2>&1

    # Convert the output sam file to bam file, sort and index the bam file
    samtools view -S -b -@ ${task.cpus} -F 4 ${base}_${ref_id}_${ref_tag}_map_ref.sam | samtools sort -@ ${task.cpus} - > ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam
    samtools index -@ ${task.cpus} ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam
    rm ${base}_${ref_id}_${ref_tag}_map_ref.sam

    if [[ ${params.dedup} == true ]]
    then
    java -Xmx${task.memory.giga}g -jar \$PICARD MarkDuplicates -I ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam -O ${base}_${ref_id}_${ref_tag}_map_ref_deduplicated.bam -M ${base}_${ref_id}_${ref_tag}_picard_output.txt --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT --REMOVE_DUPLICATES true
    # remove pre-deduplicated bam file and rename deduplicated bam file
    mv ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam ${base}_${ref_id}_${ref_tag}_map_ref_og.sorted.bam
    samtools sort -@ ${task.cpus} -T ${base}_${ref_id}_${ref_tag}_map_ref_deduplicated -o ${base}_${ref_id}_${ref_tag}_map_ref_deduplicated.sorted.bam ${base}_${ref_id}_${ref_tag}_map_ref_deduplicated.bam
    mv ${base}_${ref_id}_${ref_tag}_map_ref_deduplicated.sorted.bam ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam
    samtools index -@ ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam
    # convert deduplicated reads in bam to fastq
    samtools fastq -@ ${task.cpus} -n -0 ${base}.trimmed.fastq.gz ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam
    else
    # Convert bam of only mapped read (-F 4) to fastq
    samtools fastq -@ ${task.cpus} -n -0 ${base}.trimmed.fastq.gz ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam
    fi

    # Calling Consensus
    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}_${ref_id}_${ref_tag}.fa \\
        --min-BQ 20 \\
        --output ${base}_${ref_id}_${ref_tag}_1.mpileup \\
        ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam
    cat ${base}_${ref_id}_${ref_tag}_1.mpileup | ivar consensus -q ${params.q} -t ${params.t} -m ${params.d} -n N -p ${base}_${ref_id}_${ref_tag}.consensus1 -i ${base}_${ref_id}_${ref_tag}.consensus1

    # Get rid of leading and trailing repeated Ns and ns
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}_${ref_id}_${ref_tag}.consensus1.fa > ${base}_${ref_id}_${ref_tag}.consensus1.temp.fa
    cp ${base}_${ref_id}_${ref_tag}.consensus1.temp.fa ${base}_${ref_id}_${ref_tag}.consensus1.fa
    rm ${base}_${ref_id}_${ref_tag}.consensus1.temp.fa
 
    # Map reads to consensus1 and create bam and sorted bam files
    bbmap.sh \\
	in=${base}.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_tag}_map1.sam \\
	ref=${base}_${ref_id}_${ref_tag}.consensus1.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=80 ambiguous=random -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_tag}_map1_stats.txt 2>&1
    samtools view -S -b -@ ${task.cpus} -F 4 ${base}_${ref_id}_${ref_tag}_map1.sam | samtools sort -@ ${task.cpus} - > ${base}_${ref_id}_${ref_tag}_map1.sorted.bam
    rm ${base}_${ref_id}_${ref_tag}_map1.sam

    # second consensus generation 
    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}_${ref_id}_${ref_tag}.consensus1.fa \\
        --min-BQ 20 \\
        --output ${base}_${ref_id}_${ref_tag}_2.mpileup \\
        ${base}_${ref_id}_${ref_tag}_map1.sorted.bam
    cat ${base}_${ref_id}_${ref_tag}_2.mpileup | ivar consensus -q ${params.q} -t ${params.t} -m ${params.d} -n N -p ${base}_${ref_id}_${ref_tag}.consensus2 -i ${base}_${ref_id}_${ref_tag}.consensus2

    # Get rid of repeated Ns and ns
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}_${ref_id}_${ref_tag}.consensus2.fa > ${base}_${ref_id}_${ref_tag}.consensus2.temp.fa
    cp ${base}_${ref_id}_${ref_tag}.consensus2.temp.fa ${base}_${ref_id}_${ref_tag}.consensus2.fa
    rm ${base}_${ref_id}_${ref_tag}.consensus2.temp.fa
    
    # Align reads to consensus2 and create bam and sorted bam files 
    bbmap.sh \\
	in=${base}.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_tag}_map2.sam \\
	ref=${base}_${ref_id}_${ref_tag}.consensus2.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=80 ambiguous=random -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_tag}_map2_stats.txt 2>&1
    samtools view -S -b -@ ${task.cpus} -F 4 ${base}_${ref_id}_${ref_tag}_map2.sam | samtools sort -@ ${task.cpus} - > ${base}_${ref_id}_${ref_tag}_map2.sorted.bam
    rm ${base}_${ref_id}_${ref_tag}_map2.sam

    # Final Consensus Generation
    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}_${ref_id}_${ref_tag}.consensus2.fa \\
        --min-BQ 20 \\
        --output ${base}_${ref_id}_${ref_tag}_final.mpileup \\
        ${base}_${ref_id}_${ref_tag}_map2.sorted.bam
    cat ${base}_${ref_id}_${ref_tag}_final.mpileup | ivar consensus -q ${params.q} -t ${params.t} -m ${params.d} -n N -p ${base}_${ref_id}_${ref_tag}.consensus_final -i ${base}_${ref_id}_${ref_tag}.consensus_final

    # Get rid of repeated Ns and ns
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}_${ref_id}_${ref_tag}.consensus_final.fa > ${base}_${ref_id}_${ref_tag}.consensus_final.temp.fa
    cp ${base}_${ref_id}_${ref_tag}.consensus_final.temp.fa ${base}_${ref_id}_${ref_tag}.consensus_final.fa
    rm ${base}_${ref_id}_${ref_tag}.consensus_final.temp.fa
 
    # Align reads to final consensus and create bam and sorted bam files 
    bbmap.sh \\
	in=${base}.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_tag}_mapf.sam \\
	ref=${base}_${ref_id}_${ref_tag}.consensus_final.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=80 ambiguous=random -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_tag}_mapf_stats.txt 2>&1
    samtools view -S -b -@ ${task.cpus} -F 4 ${base}_${ref_id}_${ref_tag}_mapf.sam | samtools sort -@ ${task.cpus} - > ${base}_${ref_id}_${ref_tag}_mapf.sorted.bam
    samtools index -@ ${task.cpus} ${base}_${ref_id}_${ref_tag}_mapf.sorted.bam
    rm ${base}_${ref_id}_${ref_tag}_mapf.sam

    """
}

process VCF_Generation {
    errorStrategy 'retry'
    maxRetries 1

    input:
        tuple val(base), val(ref_id), val(ref_tag), file("${base}_${ref_id}_${ref_tag}.fa"), file("${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam")

    output:
        tuple val(base), val(ref_id), val(ref_tag), file("*.gz"), file("*.tbi"), file("*stats.txt")

    publishDir "${params.outdir}/bcftools_vcf", mode: 'copy', pattern:'*.vcf.gz'
    publishDir "${params.outdir}/bcftools_vcf", mode: 'copy', pattern:'*.tbi'        
    publishDir "${params.outdir}/bcftools_vcf", mode: 'copy', pattern:'*stats.txt'        

    script:

    """
    echo "${base}" > sample_name.list

    bcftools mpileup \\
        --fasta-ref ${base}_${ref_id}_${ref_tag}.fa \\
        --ignore-overlaps \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 0 \\
        --min-BQ 20 \\
        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
        ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam \\
        | bcftools call \\
            --output-type v \\
            --ploidy 1 \\
            --keep-alts \\
            --keep-masked-ref \\
            --multiallelic-caller \\
            --variants-only \\
        | bcftools reheader \\
            --samples sample_name.list \\
        | bcftools view \\
            --output-file ${base}.vcf.gz \\
            --output-type z \\
            --include 'INFO/DP>=10'

    tabix -p vcf -f ${base}.vcf.gz

    bcftools stats ${base}.vcf.gz > ${base}.bcftools_stats.txt

    """

}
process Serotyping {
    container 'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0'
    errorStrategy 'retry'
    maxRetries 1

    input:
	tuple val(base), val(ref_id), val(ref_tag), file("${base}_${ref_id}_${ref_tag}.consensus_final.fa")
        file BLASTDB_ALL_1
        file BLASTDB_ALL_2
        file BLASTDB_ALL_3
        file BLASTDB_ALL_4
        file BLASTDB_ALL_5
        file BLASTDB_ALL_6
        file BLASTDB_ALL_7
        file BLASTDB_ALL_8

    output:
        tuple val(base), val(ref_id), val(ref_tag)
        file("*.txt")

    publishDir "${params.outdir}/serotype", mode: 'copy', pattern:'*.txt'

    script:

    """
    #!/bin/bash

    blastn -out ${base}_${ref_id}_${ref_tag}_blast_output.txt -query ${base}_${ref_id}_${ref_tag}.consensus_final.fa -db ${BLASTDB_ALL_1} -outfmt 6 -task blastn -max_target_seqs 1 -evalue 1e-5
    # outfmt6 default values: 'qaccver saccver pident length mismatch gapopen qstart qend sstart send' 

    serotype=\$(awk 'FNR==1{print val,\$2}' ${base}_${ref_id}_${ref_tag}_blast_output.txt | cut -d "_" -f2- | cut -d "/" -f2-)
    echo \$serotype > ${base}_${ref_id}_${ref_tag}_serotype.txt

    """
}

process Summary_Generation {
    errorStrategy 'retry'
    maxRetries 1
    
    input:
        tuple val(base), val(ref_id), val(ref_tag)

    output:
	file("${base}_${ref_id}_${ref_tag}_summary.tsv")

    script:
    """
    #!/bin/bash
 
    # find real absolute path of outdir
    process_work_dir=\$PWD
    cd $workflow.launchDir
    outdir_realpath=\$(realpath ${params.outdir})
    cd \$process_work_dir

    # summary header
    echo "sample name\traw reads/pairs\tsurviving reads/pairs\treference accession\treference tag\treference header\treference length\treference num Ns\t%ref coverage\tmedian coverage\tconsensus length\tmapped reads\t%reads on target\tnum Ns\t%N\tserotype" > ${base}_${ref_id}_${ref_tag}_summary.tsv

    # get the number of total reads/pairs and suviving reads/pairs
    num_untrimmed=\$(cat \${outdir_realpath}/trim_stats/${base}_trim_stats.txt | grep "Input Read" | cut -d ":" -f2 | awk '{print \$1}')
    num_trimmed_pct=\$(cat \${outdir_realpath}/trim_stats/${base}_trim_stats.txt | grep "Input Read" | cut -d ":" -f3 | awk '{print \$1, \$2}')
    num_trimmed=\$(echo \$num_trimmed_pct | awk '{print \$1}')

    # get the reference header info
    ref_tag=\$(cat \${outdir_realpath}/viral_identification/${base}_${ref_id}_${ref_tag}_vid.txt | cut -f3)
    ref_header=\$(cat \${outdir_realpath}/viral_identification/${base}_${ref_id}_${ref_tag}_vid.txt | cut -f4)

    # get the reference genome (map_all) id, size, coverage
    ref_length=\$(grep ${ref_id} \${outdir_realpath}/bbmap_covstats_all_ref/${base}_map_all_bbmap_covstats.tsv | cut -f3)
    ref_coverage=\$(grep ${ref_id} \${outdir_realpath}/bbmap_covstats_all_ref/${base}_map_all_bbmap_covstats.tsv | cut -f5)
    median_coverage=\$(grep ${ref_id} \${outdir_realpath}/bbmap_covstats_all_ref/${base}_map_all_bbmap_covstats.tsv | cut -f10)

    # get the consensus final length
    consensus_length=\$(awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length(\$0)}END{print l}' \${outdir_realpath}/consensus_final/${base}_${ref_id}_${ref_tag}.consensus_final.fa | awk 'FNR==2{print val,\$1}')

    # get the number and the percentage of trimmed reads mapped to the reference genome (map_ref)
    mapped_reads=\$(cat \${outdir_realpath}/map_consensus_final_stats/${base}_${ref_id}_${ref_tag}_mapf_stats.txt | grep "mapped:" | cut -f3)
    # works for single-end and paired-end reads
    mapped_reads=\$(echo \$mapped_reads | awk '{print \$1+\$2}')
 
    # calculate % reads on target for single-end and paired-end
    if [[ ${params.pe} == false ]]
    then
    percent_mapped_reads=\$(echo "\$mapped_reads/\$num_trimmed*100" | bc -l | awk 'FNR==1{print val,\$1}')
    else
    percent_mapped_reads=\$(echo "\$mapped_reads/\$num_trimmed/2*100" | bc -l | awk 'FNR==1{print val,\$1}')
    fi

    # get the number of Ns in the consensus final 
    num_ns=\$(grep -v "^>" \${outdir_realpath}/consensus_final/${base}_${ref_id}_${ref_tag}.consensus_final.fa | tr -cd N | wc -c | awk 'FNR==1{print val,\$1}')

    # get the number of Ns in the consensus final 
    ref_num_ns=\$(grep -v "^>" \${outdir_realpath}/ref_genome/${base}_${ref_id}_${ref_tag}.fa | tr -cd N | wc -c | awk 'FNR==1{print val,\$1}')

    # get the percentage of Ns in the consensus final
    percent_n=\$(echo "\$num_ns/\$consensus_length*100" | bc -l | awk 'FNR==1{print val,\$1}')

    # get the serotype
    serotype=\$(awk 'FNR==1{print \$1}' \${outdir_realpath}/serotype/${base}_${ref_id}_${ref_tag}_serotype.txt)

    echo "${base}\t\$num_untrimmed\t\$num_trimmed_pct\t${ref_id}\t\$ref_tag\t\$ref_header\t\$ref_length\t\$ref_num_ns\t\$ref_coverage\t\$median_coverage\t\$consensus_length\t\$mapped_reads\t\$percent_mapped_reads\t\$num_ns\t\$percent_n\t\$serotype" >> ${base}_${ref_id}_${ref_tag}_summary.tsv

    """
}

process Final_Processing {
    errorStrategy 'retry'
    maxRetries 1

    input:
	file("*_summary.tsv")

    output: 
	file("run_summary.tsv")
	file("failed_assembly_summary.tsv") optional true


    publishDir "${params.outdir}", mode: 'copy', pattern:'run_summary.tsv'
    publishDir "${params.outdir}", mode: 'copy', pattern:'failed_assembly_summary.tsv'

    script:
    """
    #!/bin/bash

    # generate a summary on samples that pass median coverage threshold and have a consensus genome generated.
    echo "sample name\traw reads/pairs\tsurviving reads/pairs\treference accession\treference tag\treference header\treference length\treference num Ns\t%ref coverage\tmedian coverage\tconsensus length\tmapped reads\t%reads on target\tnum Ns\t%N\tserotype" > summary.tsv
    awk '(NR == 2) || (FNR > 1)' *_summary.tsv >> summary.tsv 
    head -1 summary.tsv > run_summary.tsv
    awk 'NR>1' < summary.tsv | sort -k1 >> run_summary.tsv

    # find real absolute path of outdir
    process_work_dir=\$PWD
    cd $workflow.launchDir
    outdir_realpath=\$(realpath ${params.outdir})
    cd \$process_work_dir

    # summary stats for failed assembly
    if [ -n "\$(ls -A \${outdir_realpath}/failed_assembly 2>/dev/null)" ]
    then
	cp \${outdir_realpath}/failed_assembly/*_failed_assembly.tsv .
	cat *_failed_assembly.tsv | awk 'NR==1 || NR % 2 == 0' > failed_assembly_summary.tsv
    fi
	
    """
}

process Trimming_PE { 
    errorStrategy 'ignore'

    input:
        tuple val(base), file(R1), file(R2) // from input_read_ch
        file ADAPTERS

    output: 
        tuple val(base), file("${base}_R1.trimmed.fastq.gz"), file("${base}_R2.trimmed.fastq.gz")
	file("*_trim_stats.txt")

    publishDir "${params.outdir}/trimmed_fastq", mode: 'copy', pattern:'*.trimmed.fastq.gz'
    publishDir "${params.outdir}/trim_stats", mode: 'copy', pattern:'*_trim_stats.txt'

    script:
    """
    #!/bin/bash

    if [[ ${R1} == *.fastq && ${R2} == *.fastq ]]
    then
	if [[ ${params.sample} == false ]]
	then
            gzip ${R1} && gzip ${R2}
    	    java -Xmx${task.memory.giga}g -jar \$TRIMMOMATIC PE -threads ${task.cpus} ${R1}.gz ${R2}.gz ${base}_R1.trimmed.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}_R2.trimmed.fastq.gz ${base}.R2.unpaired.fastq.gz \
    	    ILLUMINACLIP:${ADAPTERS}:${params.SETTING} LEADING:${params.LEADING} TRAILING:${params.TRAILING} SLIDINGWINDOW:${params.SWINDOW} MINLEN:${params.MINLEN} &> ${base}_trim_stats.txt
	else
	    seqtk sample -s 100 ${R1} ${params.sample} | gzip > ${base}_R1_sampled.fastq.gz
	    seqtk sample -s 100 ${R2} ${params.sample} | gzip > ${base}_R2_sampled.fastq.gz
    	    java -Xmx${task.memory.giga}g -jar \$TRIMMOMATIC PE -threads ${task.cpus} ${base}_R1_sampled.fastq.gz ${base}_R2_sampled.fastq.gz ${base}_R1.trimmed.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}_R2.trimmed.fastq.gz ${base}.R2.unpaired.fastq.gz \
    	    ILLUMINACLIP:${ADAPTERS}:${params.SETTING} LEADING:${params.LEADING} TRAILING:${params.TRAILING} SLIDINGWINDOW:${params.SWINDOW} MINLEN:${params.MINLEN} &> ${base}_trim_stats.txt
	fi
    elif [[ ${R1} == *.fastq.gz && ${R2} == *.fastq.gz ]]
    then 
	if [[ ${params.sample} == false ]]
	then
    	    java -Xmx${task.memory.giga}g -jar \$TRIMMOMATIC PE -threads ${task.cpus} ${R1} ${R2} ${base}_R1.trimmed.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}_R2.trimmed.fastq.gz ${base}.R2.unpaired.fastq.gz \
    	    ILLUMINACLIP:${ADAPTERS}:${params.SETTING} LEADING:${params.LEADING} TRAILING:${params.TRAILING} SLIDINGWINDOW:${params.SWINDOW} MINLEN:${params.MINLEN} &> ${base}_trim_stats.txt
	else
	    seqtk sample -s 100 ${R1} ${params.sample} | gzip > ${base}_R1_sampled.fastq.gz
	    seqtk sample -s 100 ${R2} ${params.sample} | gzip > ${base}_R2_sampled.fastq.gz
    	    java -Xmx${task.memory.giga}g -jar \$TRIMMOMATIC PE -threads ${task.cpus} ${base}_R1_sampled.fastq.gz ${base}_R2_sampled.fastq.gz ${base}_R1.trimmed.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}_R2.trimmed.fastq.gz ${base}.R2.unpaired.fastq.gz \
    	    ILLUMINACLIP:${ADAPTERS}:${params.SETTING} LEADING:${params.LEADING} TRAILING:${params.TRAILING} SLIDINGWINDOW:${params.SWINDOW} MINLEN:${params.MINLEN} &> ${base}_trim_stats.txt
	fi
    fi

    """
}

process Aligning_PE {
    errorStrategy 'retry'
    maxRetries 1

    input: 
        tuple val(base), file("${base}_R1.trimmed.fastq.gz"), file("${base}_R2.trimmed.fastq.gz")
        file ref

    output:
	tuple val(base), file("${base}_map_all_bbmap_covstats.tsv"), file("${base}_map_all_stats.txt")

    publishDir "${params.outdir}/bbmap_covstats_all_ref", mode: 'copy', pattern:'*_map_all_bbmap_covstats.tsv'
    
    script:

    """
    #!/bin/bash

    # Map to the multifasta reference
    bbmap.sh \\
        in=${base}_R1.trimmed.fastq.gz \\
	in2=${base}_R2.trimmed.fastq.gz \\
        outm=${base}_map_all.sam \\
        ref=${ref} \\
        threads=${task.cpus} \\
        covstats=${base}_map_all_bbmap_covstats.tsv \\
        local=true interleaved=false maxindel=80 -Xmx${task.memory.giga}g > ${base}_map_all_stats.txt 2>&1

    # sort bbmap_out.txt based on median_fold (column 10) in descending order 
    head -1 ${base}_map_all_bbmap_covstats.tsv > ${base}_map_all_bbmap_covstats.temp.txt
    awk 'NR>1' < ${base}_map_all_bbmap_covstats.tsv | sort -t \$'\t' -nrk10 >> ${base}_map_all_bbmap_covstats.temp.txt
    mv ${base}_map_all_bbmap_covstats.temp.txt ${base}_map_all_bbmap_covstats.tsv

    """
}

process Consensus_Generation_Prep_PE {
    errorStrategy 'retry'
    maxRetries 1

    input:
	file "*_vid.txt"
        file ref

    output:
	tuple env(base), env(ref_id), env(ref_tag), file("*.fa") 

    publishDir "${params.outdir}/ref_genome", mode: 'copy', pattern: '*.fa'

    script:
    """
    #!/bin/bash

    base=\$(head -1 *_vid.txt | cut -f1)
    ref_id=\$(head -1 *_vid.txt | cut -f2)
    ref_tag=\$(head -1 *_vid.txt | cut -f3)
    samtools faidx ${ref} \$ref_id > \${base}'_'\${ref_id}'_'\${ref_tag}'.fa'     

    """
}

process Consensus_Generation_PE {
    errorStrategy 'retry'
    maxRetries 1 

    input:
	tuple val(base), val(ref_id), val(ref_tag), file("${base}_${ref_id}_${ref_tag}.fa") 
    output:
	tuple val(base), val(ref_id), val(ref_tag), file("${base}_${ref_id}_${ref_tag}.fa"), file("${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam")
        tuple val(base), val(ref_id), val(ref_tag), file("${base}_${ref_id}_${ref_tag}.consensus_final.fa")
	tuple file("*.txt"), file("*.fa"), file("*bam*"), file("*.bai") 

    publishDir "${params.outdir}/map_ref_bam_sorted", mode: 'copy', pattern: '*_map_ref.sorted.bam'
    publishDir "${params.outdir}/map_ref_bam_sorted", mode: 'copy', pattern: '*_map_ref.sorted.bam.bai'
    publishDir "${params.outdir}/consensus_final", mode: 'copy', pattern:'*.consensus_final.fa' 
    publishDir "${params.outdir}/consensus_final_bam_sorted", mode: 'copy', pattern:'*_mapf.sorted.bam' 
    publishDir "${params.outdir}/consensus_final_bam_sorted", mode: 'copy', pattern:'*_mapf.sorted.bam.bai' 
    publishDir "${params.outdir}/map_consensus_final_stats", mode: 'copy', pattern: '*_mapf_stats.txt'

    script:
    """
    #!/bin/bash
   
    # find real absolute path of outdir
    process_work_dir=\$PWD
    cd $workflow.launchDir
    outdir_realpath=\$(realpath ${params.outdir})
    cd \$process_work_dir

    # Map the reads to the reference
    bbmap.sh \\
	in=\${outdir_realpath}/trimmed_fastq/${base}_R1.trimmed.fastq.gz \\
	in2=\${outdir_realpath}/trimmed_fastq/${base}_R2.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_tag}_map_ref.sam \\
	ref=${base}_${ref_id}_${ref_tag}.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=80 ambiguous=random -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_tag}_map_ref_stats.txt 2>&1

    # Convert the output sam file to bam file, sort and index the bam file
    samtools view -S -b -@ ${task.cpus} -F 4 ${base}_${ref_id}_${ref_tag}_map_ref.sam | samtools sort -@ ${task.cpus} - > ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam
    samtools index -@ ${task.cpus} ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam
    rm ${base}_${ref_id}_${ref_tag}_map_ref.sam

    if [[ ${params.dedup} == true ]]
    then
    java -Xmx${task.memory.giga}g -jar \$PICARD MarkDuplicates -I ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam -O ${base}_${ref_id}_${ref_tag}_map_ref_deduplicated.bam -M ${base}_${ref_id}_${ref_tag}_picard_output.txt --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT --REMOVE_DUPLICATES true
    mv ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam ${base}_${ref_id}_${ref_tag}_map_ref_og.sorted.bam
    samtools sort -@ ${task.cpus} -T ${base}_${ref_id}_${ref_tag}_map_ref_deduplicated -o ${base}_${ref_id}_${ref_tag}_map_ref_deduplicated.sorted.bam ${base}_${ref_id}_${ref_tag}_map_ref_deduplicated.bam
    mv ${base}_${ref_id}_${ref_tag}_map_ref_deduplicated.sorted.bam ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam
    samtools index -@ ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam
    # convert deduplicated reads in bam to fastq
    samtools collate -@ ${task.cpus} -O ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam | samtools fastq -@ ${task.cpus} -1 ${base}_R1.trimmed.fastq.gz -2 ${base}_R2.trimmed.fastq.gz -0 /dev/null -s /dev/null -n 
    else
    # convert mapped reads from bam to fastq
    samtools collate -@ ${task.cpus} -O ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam | samtools fastq -@ ${task.cpus} -1 ${base}_R1.trimmed.fastq.gz -2 ${base}_R2.trimmed.fastq.gz -0 /dev/null -s /dev/null -n 
    fi

    # Calling Consensus
    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}_${ref_id}_${ref_tag}.fa \\
        --min-BQ 20 \\
        --output ${base}_${ref_id}_${ref_tag}_1.mpileup \\
        ${base}_${ref_id}_${ref_tag}_map_ref.sorted.bam
    cat ${base}_${ref_id}_${ref_tag}_1.mpileup | ivar consensus -q ${params.q} -t ${params.t} -m ${params.d} -n N -p ${base}_${ref_id}_${ref_tag}.consensus1 -i ${base}_${ref_id}_${ref_tag}.consensus1

    # Get rid of leading and trailing repeated Ns and ns
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}_${ref_id}_${ref_tag}.consensus1.fa > ${base}_${ref_id}_${ref_tag}.consensus1.temp.fa
    mv ${base}_${ref_id}_${ref_tag}.consensus1.temp.fa ${base}_${ref_id}_${ref_tag}.consensus1.fa

    # Map reads to consensus1 and create bam and sorted bam files
    bbmap.sh \\
	in=${base}_R1.trimmed.fastq.gz \\
	in2=${base}_R2.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_tag}_map1.sam \\
	ref=${base}_${ref_id}_${ref_tag}.consensus1.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=80 ambiguous=random -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_tag}_map1_stats.txt 2>&1
    samtools view -S -b -@ ${task.cpus} -F 4 ${base}_${ref_id}_${ref_tag}_map1.sam | samtools sort -@ ${task.cpus} - > ${base}_${ref_id}_${ref_tag}_map1.sorted.bam
    rm ${base}_${ref_id}_${ref_tag}_map1.sam

    # second consensus generation 
    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}_${ref_id}_${ref_tag}.consensus1.fa \\
        --min-BQ 20 \\
        --output ${base}_${ref_id}_${ref_tag}_2.mpileup \\
        ${base}_${ref_id}_${ref_tag}_map1.sorted.bam
    cat ${base}_${ref_id}_${ref_tag}_2.mpileup | ivar consensus -q ${params.q} -t ${params.t} -m ${params.d} -n N -p ${base}_${ref_id}_${ref_tag}.consensus2 -i ${base}_${ref_id}_${ref_tag}.consensus2

    # Get rid of repeated Ns and ns
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}_${ref_id}_${ref_tag}.consensus2.fa > ${base}_${ref_id}_${ref_tag}.consensus2.temp.fa
    mv ${base}_${ref_id}_${ref_tag}.consensus2.temp.fa ${base}_${ref_id}_${ref_tag}.consensus2.fa

    # Align reads to consensus2 and create bam and sorted bam files 
    bbmap.sh \\
	in=${base}_R1.trimmed.fastq.gz \\
	in2=${base}_R2.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_tag}_map2.sam \\
	ref=${base}_${ref_id}_${ref_tag}.consensus2.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=80 ambiguous=random -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_tag}_map2_stats.txt 2>&1
    samtools view -S -b -@ ${task.cpus} -F 4 ${base}_${ref_id}_${ref_tag}_map2.sam | samtools sort -@ {$task.cpus} - > ${base}_${ref_id}_${ref_tag}_map2.sorted.bam
    rm ${base}_${ref_id}_${ref_tag}_map2.sam

    # Final Consensus Generation
    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}_${ref_id}_${ref_tag}.consensus2.fa \\
        --min-BQ 20 \\
        --output ${base}_${ref_id}_${ref_tag}_final.mpileup \\
        ${base}_${ref_id}_${ref_tag}_map2.sorted.bam
    cat ${base}_${ref_id}_${ref_tag}_final.mpileup | ivar consensus -q ${params.q} -t ${params.t} -m ${params.d} -n N -p ${base}_${ref_id}_${ref_tag}.consensus_final -i ${base}_${ref_id}_${ref_tag}.consensus_final

    # Get rid of repeated Ns and ns
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}_${ref_id}_${ref_tag}.consensus_final.fa > ${base}_${ref_id}_${ref_tag}.consensus_final.temp.fa
    mv ${base}_${ref_id}_${ref_tag}.consensus_final.temp.fa ${base}_${ref_id}_${ref_tag}.consensus_final.fa

    # Align reads to final consensus and create bam and sorted bam files 
    bbmap.sh \\
	in=${base}_R1.trimmed.fastq.gz \\
	in2=${base}_R2.trimmed.fastq.gz \\
	outm=${base}_${ref_id}_${ref_tag}_mapf.sam \\
	ref=${base}_${ref_id}_${ref_tag}.consensus_final.fa \\
	threads=${task.cpus} \\
	local=true interleaved=false maxindel=80 ambiguous=random -Xmx${task.memory.giga}g > ${base}_${ref_id}_${ref_tag}_mapf_stats.txt 2>&1
    samtools view -S -b -@ ${task.cpus} -F 4 ${base}_${ref_id}_${ref_tag}_mapf.sam | samtools sort -@ ${task.cpus} - > ${base}_${ref_id}_${ref_tag}_mapf.sorted.bam
    samtools index -@ ${task.cpus} ${base}_${ref_id}_${ref_tag}_mapf.sorted.bam
    rm ${base}_${ref_id}_${ref_tag}_mapf.sam
    """
}
