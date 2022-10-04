#!/bin/bash

   # find real absolute path of outdir
   process_work_dir=$PWD
   cd /Users/greningerlab/AS/revica
   outdir_realpath=$(realpath test)
   cd $process_work_dir

   # Map the reads to the reference
   bbmap.sh \
in=${outdir_realpath}/trimmed_fastq/example.trimmed.fastq.gz \
outm=example_KF923898.1_HCOV_map_ref.sam \
ref=example_KF923898.1_HCOV.fa \
threads=8 \
local=true interleaved=false maxindel=80 ambiguous=random -Xmx8g > example_KF923898.1_HCOV_map_ref_stats.txt 2>&1

   # Convert the output sam file to bam file, sort and index the bam file
   samtools view -S -b -@ 8 -F 4 example_KF923898.1_HCOV_map_ref.sam | samtools sort -@ 8 - > example_KF923898.1_HCOV_map_ref.sorted.bam
   samtools index -@ 8 example_KF923898.1_HCOV_map_ref.sorted.bam
   rm example_KF923898.1_HCOV_map_ref.sam

   if [[ false == true ]]
   then
   java -Xmx8g -jar $PICARD MarkDuplicates -I example_KF923898.1_HCOV_map_ref.sorted.bam -O example_KF923898.1_HCOV_map_ref_deduplicated.bam -M example_KF923898.1_HCOV_picard_output.txt --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT --REMOVE_DUPLICATES true
   # remove pre-deduplicated bam file and rename deduplicated bam file
   mv example_KF923898.1_HCOV_map_ref.sorted.bam example_KF923898.1_HCOV_map_ref_og.sorted.bam
   samtools sort -@ 8 -T example_map_ref_deduplicated -o example_map_ref_deduplicated.sorted.bam example_map_ref_deduplicated.bam
   mv example_KF923898.1_HCOV_map_ref_deduplicated.sorted.bam example_KF923898.1_HCOV_map_ref.sorted.bam
   samtools index -@ example_KF923898.1_HCOV_map_ref.sorted.bam
   # convert deduplicated reads in bam to fastq
   samtools fastq -@ 8 -n -0 example.trimmed.fastq.gz example_KF923898.1_HCOV_map_ref.sorted.bam
   else
   # Convert bam of only mapped read (-F 4) to fastq
   samtools fastq -@ 8 -n -0 example.trimmed.fastq.gz example_KF923898.1_HCOV_map_ref.sorted.bam
   fi

   # Calling Consensus
   samtools mpileup \
       --count-orphans \
       --no-BAQ \
       --max-depth 50000 \
       --fasta-ref example_KF923898.1_HCOV.fa \
       --min-BQ 20 \
       --output example_KF923898.1_HCOV_1.mpileup \
       example_KF923898.1_HCOV_map_ref.sorted.bam
   cat example_KF923898.1_HCOV_1.mpileup | ivar consensus -q 20 -t 0.6 -m 5 -n N -p example_KF923898.1_HCOV.consensus1 -i example_KF923898.1_HCOV.consensus1

   # Get rid of leading and trailing repeated Ns and ns
   seqkit -is replace -p "^n+|n+$" -r "" example_KF923898.1_HCOV.consensus1.fa > example_KF923898.1_HCOV.consensus1.temp.fa
   cp example_KF923898.1_HCOV.consensus1.temp.fa example_KF923898.1_HCOV.consensus1.fa
   rm example_KF923898.1_HCOV.consensus1.temp.fa

   # Map reads to consensus1 and create bam and sorted bam files
   bbmap.sh \
in=example.trimmed.fastq.gz \
outm=example_KF923898.1_HCOV_map1.sam \
ref=example_KF923898.1_HCOV.consensus1.fa \
threads=8 \
local=true interleaved=false maxindel=80 ambiguous=random -Xmx8g > example_KF923898.1_HCOV_map1_stats.txt 2>&1
   samtools view -S -b -@ 8 -F 4 example_KF923898.1_HCOV_map1.sam | samtools sort -@ 8 - > example_KF923898.1_HCOV_map1.sorted.bam
   rm example_KF923898.1_HCOV_map1.sam

   # second consensus generation 
   samtools mpileup \
       --count-orphans \
       --no-BAQ \
       --max-depth 50000 \
       --fasta-ref example_KF923898.1_HCOV.consensus1.fa \
       --min-BQ 20 \
       --output example_KF923898.1_HCOV_2.mpileup \
       example_KF923898.1_HCOV_map1.sorted.bam
   cat example_KF923898.1_HCOV_2.mpileup | ivar consensus -q 20 -t 0.6 -m 5 -n N -p example_KF923898.1_HCOV.consensus2 -i example_KF923898.1_HCOV.consensus2

   # Get rid of repeated Ns and ns
   seqkit -is replace -p "^n+|n+$" -r "" example_KF923898.1_HCOV.consensus2.fa > example_KF923898.1_HCOV.consensus2.temp.fa
   cp example_KF923898.1_HCOV.consensus2.temp.fa example_KF923898.1_HCOV.consensus2.fa
   rm example_KF923898.1_HCOV.consensus2.temp.fa

   # Align reads to consensus2 and create bam and sorted bam files 
   bbmap.sh \
in=example.trimmed.fastq.gz \
outm=example_KF923898.1_HCOV_map2.sam \
ref=example_KF923898.1_HCOV.consensus2.fa \
threads=8 \
local=true interleaved=false maxindel=80 ambiguous=random -Xmx8g > example_KF923898.1_HCOV_map2_stats.txt 2>&1
   samtools view -S -b -@ 8 -F 4 example_KF923898.1_HCOV_map2.sam | samtools sort -@ 8 - > example_KF923898.1_HCOV_map2.sorted.bam
   rm example_KF923898.1_HCOV_map2.sam

   # Final Consensus Generation
   samtools mpileup \
       --count-orphans \
       --no-BAQ \
       --max-depth 50000 \
       --fasta-ref example_KF923898.1_HCOV.consensus2.fa \
       --min-BQ 20 \
       --output example_KF923898.1_HCOV_final.mpileup \
       example_KF923898.1_HCOV_map2.sorted.bam
   cat example_KF923898.1_HCOV_final.mpileup | ivar consensus -q 20 -t 0.6 -m 5 -n N -p example_KF923898.1_HCOV.consensus_final -i example_KF923898.1_HCOV.consensus_final

   # Get rid of repeated Ns and ns
   seqkit -is replace -p "^n+|n+$" -r "" example_KF923898.1_HCOV.consensus_final.fa > example_KF923898.1_HCOV.consensus_final.temp.fa
   cp example_KF923898.1_HCOV.consensus_final.temp.fa example_KF923898.1_HCOV.consensus_final.fa
   rm example_KF923898.1_HCOV.consensus_final.temp.fa

   # Align reads to final consensus and create bam and sorted bam files 
   bbmap.sh \
in=example.trimmed.fastq.gz \
outm=example_KF923898.1_HCOV_mapf.sam \
ref=example_KF923898.1_HCOV.consensus_final.fa \
threads=8 \
local=true interleaved=false maxindel=80 ambiguous=random -Xmx8g > example_KF923898.1_HCOV_mapf_stats.txt 2>&1
   samtools view -S -b -@ 8 -F 4 example_KF923898.1_HCOV_mapf.sam | samtools sort -@ 8 - > example_KF923898.1_HCOV_mapf.sorted.bam
   samtools index -@ 8 example_KF923898.1_HCOV_mapf.sorted.bam
   rm example_KF923898.1_HCOV_mapf.sam
