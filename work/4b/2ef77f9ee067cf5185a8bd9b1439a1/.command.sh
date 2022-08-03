#!/bin/bash

   # generate a summary on samples that pass median coverage threshold and have a consensus genome generated.
   echo "sample name	raw reads/pairs	surviving reads/pairs	reference accession	reference tag	reference header	reference length	reference num Ns	%ref coverage	median coverage	consensus length	mapped reads	%reads on target	num Ns	%N	serotype" > summary.tsv
   awk '(NR == 2) || (FNR > 1)' *_summary.tsv >> summary.tsv 
   head -1 summary.tsv > run_summary.tsv
   awk 'NR>1' < summary.tsv | sort -k1 >> run_summary.tsv

   # find real absolute path of outdir
   process_work_dir=$PWD
   cd /Users/greningerlab/AS/revica
   outdir_realpath=$(realpath /Users/greningerlab/Downloads/revica_test)
   cd $process_work_dir

   # summary stats for failed assembly
   if [ -n "$(ls -A ${outdir_realpath}/failed_assembly 2>/dev/null)" ]
   then
cp ${outdir_realpath}/failed_assembly/*_failed_assembly.tsv .
cat *_failed_assembly.tsv | awk 'NR==1 || NR % 2 == 0' > failed_assembly_summary.tsv
   fi
