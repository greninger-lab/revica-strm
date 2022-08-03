#!/bin/bash

   if [[ example.fastq.gz == *.fastq ]]
   then
   base=$(basename example.fastq.gz ".fastq")
   	if [[ false == false ]]
   	then 
   	    trimmomatic SE -threads 2 example.fastq.gz $base.trimmed.fastq ILLUMINACLIP:adapters.fa:2:30:10:1:true     	    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35 &> ${base}_trim_stats.txt
   	    gzip $base.trimmed.fastq
else
    seqtk sample example.fastq.gz false > ${base}_sampled.fastq
   	    trimmomatic SE -threads 2 ${base}_sampled.fastq $base.trimmed.fastq ILLUMINACLIP:adapters.fa:2:30:10:1:true     	    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35 &> ${base}_trim_stats.txt
   	    gzip $base.trimmed.fastq
fi
   elif [[ example.fastq.gz == *.fastq.gz ]] 
   then
   base=$(basename example.fastq.gz ".fastq.gz")
if [[ false == false ]]
then
   	    trimmomatic SE -threads 2 example.fastq.gz $base.trimmed.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10:1:true     	    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35 &> ${base}_trim_stats.txt
else
    seqtk sample example.fastq.gz false > ${base}_sampled.fastq
    gzip ${base}_sampled.fastq
   	    trimmomatic SE -threads 2 ${base}_sampled.fastq.gz ${base}.trimmed.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10:1:true     	    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35 &> ${base}_trim_stats.txt
fi
   fi

# capture process environment
set +u
echo base=$base > .command.env
