#!/bin/bash

   # find real absolute path of outdir
   process_work_dir=$PWD
   cd /Users/greningerlab/AS/revica
   outdir_realpath=$(realpath /Users/greningerlab/Downloads/revica_test)
   cd $process_work_dir

   # Map the reads to the reference
   bbmap.sh \
in=${outdir_realpath}/trimmed_fastq/example.trimmed.fastq.gz \
outm=example_KY369892.1_RV_map_ref.sam \
ref=example_KY369892.1_RV.fa \
threads=8 \
local=true interleaved=false maxindel=80 ambiguous=random -Xmx8g > example_KY369892.1_RV_map_ref_stats.txt 2>&1

   # Convert the output sam file to bam file, sort and index the bam file
   samtools view -S -b -@ 8 -F 4 example_KY369892.1_RV_map_ref.sam | samtools sort -@ 8 - > example_KY369892.1_RV_map_ref.sorted.bam
   rm example_KY369892.1_RV_map_ref.sam

   if [[ false == true ]]
   then
   picard MarkDuplicates -I example_KY369892.1_RV_map_ref.sorted.bam -O example_KY369892.1_RV_map_ref_deduplicated.sorted.bam -M example_KY369892.1_RV_picard_output.txt -REMOVE_DUPLICATES true
   # remove pre-deduplicated bam file and rename deduplicated bam file
   mv example_KY369892.1_RV_map_ref.sorted.bam example_KY369892.1_RV_map_ref_og.sorted.bam
   mv example_KY369892.1_RV_map_ref_deduplicated.sorted.bam example_KY369892.1_RV_map_ref.sorted.bam
   # convert deduplicated reads in bam to fastq
   samtools fastq -@ 8 -n -0 example.trimmed.fastq.gz example_KY369892.1_RV_map_ref.sorted.bam
   else
   # Convert bam of only mapped read (-F 4) to fastq
   samtools fastq -@ 8 -n -0 example.trimmed.fastq.gz example_KY369892.1_RV_map_ref.sorted.bam
   fi

   # Calling Consensus
   samtools mpileup \
       --count-orphans \
       --no-BAQ \
       --max-depth 50000 \
       --fasta-ref example_KY369892.1_RV.fa \
       --min-BQ 15 \
       --output example_KY369892.1_RV_1.mpileup \
       example_KY369892.1_RV_map_ref.sorted.bam
   cat example_KY369892.1_RV_1.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p example_KY369892.1_RV.consensus1

   # Get rid of leading and trailing repeated Ns and ns
   seqkit -is replace -p "^n+|n+$" -r "" example_KY369892.1_RV.consensus1.fa > example_KY369892.1_RV.consensus1.temp.fa
   cp example_KY369892.1_RV.consensus1.temp.fa example_KY369892.1_RV.consensus1.fa
   rm example_KY369892.1_RV.consensus1.temp.fa

   # Modify first consensus header
   cp example_KY369892.1_RV.consensus1.fa example_KY369892.1_RV.consensus1.fa.backup
   echo '>example_KY369892.1_RV.consensus1' > example_KY369892.1_RV.consensus1.fa
   tail -n+2 example_KY369892.1_RV.consensus1.fa.backup >> example_KY369892.1_RV.consensus1.fa
   rm example_KY369892.1_RV.consensus1.fa.backup

   # Map reads to consensus1 and create bam and sorted bam files
   bbmap.sh \
in=example.trimmed.fastq.gz \
outm=example_KY369892.1_RV_map1.sam \
ref=example_KY369892.1_RV.consensus1.fa \
threads=8 \
local=true interleaved=false maxindel=80 ambiguous=random -Xmx8g > example_KY369892.1_RV_map1_stats.txt 2>&1
   samtools view -S -b -@ 8 -F 4 example_KY369892.1_RV_map1.sam | samtools sort -@ 8 - > example_KY369892.1_RV_map1.sorted.bam
   rm example_KY369892.1_RV_map1.sam

   # second consensus generation 
   samtools mpileup \
       --count-orphans \
       --no-BAQ \
       --max-depth 50000 \
       --fasta-ref example_KY369892.1_RV.consensus1.fa \
       --min-BQ 15 \
       --output example_KY369892.1_RV_2.mpileup \
       example_KY369892.1_RV_map1.sorted.bam
   cat example_KY369892.1_RV_2.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p example_KY369892.1_RV.consensus2

   # Get rid of repeated Ns and ns
   seqkit -is replace -p "^n+|n+$" -r "" example_KY369892.1_RV.consensus2.fa > example_KY369892.1_RV.consensus2.temp.fa
   cp example_KY369892.1_RV.consensus2.temp.fa example_KY369892.1_RV.consensus2.fa
   rm example_KY369892.1_RV.consensus2.temp.fa

   # Modify second consensus header
   cp example_KY369892.1_RV.consensus2.fa example_KY369892.1_RV.consensus2.fa.backup
   echo '>example_KY369892.1_RV.consensus2' > example_KY369892.1_RV.consensus2.fa
   tail -n+2 example_KY369892.1_RV.consensus2.fa.backup >> example_KY369892.1_RV.consensus2.fa
   rm example_KY369892.1_RV.consensus2.fa.backup

   # Align reads to consensus2 and create bam and sorted bam files 
   bbmap.sh \
in=example.trimmed.fastq.gz \
outm=example_KY369892.1_RV_map2.sam \
ref=example_KY369892.1_RV.consensus2.fa \
threads=8 \
local=true interleaved=false maxindel=80 ambiguous=random -Xmx8g > example_KY369892.1_RV_map2_stats.txt 2>&1
   samtools view -S -b -@ 8 -F 4 example_KY369892.1_RV_map2.sam | samtools sort -@ 8 - > example_KY369892.1_RV_map2.sorted.bam
   rm example_KY369892.1_RV_map2.sam

   # Final Consensus Generation
   samtools mpileup \
       --count-orphans \
       --no-BAQ \
       --max-depth 50000 \
       --fasta-ref example_KY369892.1_RV.consensus2.fa \
       --min-BQ 15 \
       --output example_KY369892.1_RV_final.mpileup \
       example_KY369892.1_RV_map2.sorted.bam
   cat example_KY369892.1_RV_final.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p example_KY369892.1_RV.consensus_final

   # Get rid of repeated Ns and ns
   seqkit -is replace -p "^n+|n+$" -r "" example_KY369892.1_RV.consensus_final.fa > example_KY369892.1_RV.consensus_final.temp.fa
   cp example_KY369892.1_RV.consensus_final.temp.fa example_KY369892.1_RV.consensus_final.fa
   rm example_KY369892.1_RV.consensus_final.temp.fa

   # Modify final consensus header
   cp example_KY369892.1_RV.consensus_final.fa example_KY369892.1_RV.consensus_final.fa.backup
   echo '>example_KY369892.1_RV' > example_KY369892.1_RV.consensus_final.fa
   tail -n+2 example_KY369892.1_RV.consensus_final.fa.backup >> example_KY369892.1_RV.consensus_final.fa
   rm example_KY369892.1_RV.consensus_final.fa.backup

   # Align reads to final consensus and create bam and sorted bam files 
   bbmap.sh \
in=example.trimmed.fastq.gz \
outm=example_KY369892.1_RV_mapf.sam \
ref=example_KY369892.1_RV.consensus_final.fa \
threads=8 \
local=true interleaved=false maxindel=80 ambiguous=random -Xmx8g > example_KY369892.1_RV_mapf_stats.txt 2>&1
   samtools view -S -b -@ 8 -F 4 example_KY369892.1_RV_mapf.sam | samtools sort -@ 8 - > example_KY369892.1_RV_mapf.sorted.bam
   rm example_KY369892.1_RV_mapf.sam
