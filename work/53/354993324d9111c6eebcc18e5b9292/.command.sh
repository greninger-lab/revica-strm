#!/bin/bash

blastn -out example_KF923898.1_HCOV_blast_output.txt -query example_KF923898.1_HCOV.consensus_final.fa -db VP1_164_annotated_nospaces.fasta -outfmt 6 -task blastn -max_target_seqs 1 -evalue 1e-5
# outfmt6 default values: 'qaccver saccver pident length mismatch gapopen qstart qend sstart send' 

serotype=$(awk 'FNR==1{print val,$2}' example_KF923898.1_HCOV_blast_output.txt | cut -d "_" -f2- | cut -d "/" -f2-)
echo $serotype > example_KF923898.1_HCOV_serotype.txt
