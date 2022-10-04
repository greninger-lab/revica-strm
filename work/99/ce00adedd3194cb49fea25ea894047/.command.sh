#!/bin/bash

# find real absolute path of outdir
process_work_dir=$PWD
cd /Users/greningerlab/AS/revica
outdir_realpath=$(realpath test)
cd $process_work_dir

# summary header
echo "sample name	raw reads/pairs	surviving reads/pairs	reference accession	reference tag	reference header	reference length	reference num Ns	%ref coverage	median coverage	consensus length	mapped reads	%reads on target	num Ns	%N	serotype" > example_KY369892.1_RV_summary.tsv

# get the number of total reads/pairs and suviving reads/pairs
num_untrimmed=$(cat ${outdir_realpath}/trim_stats/example_trim_stats.txt | grep "Input Read" | cut -d ":" -f2 | awk '{print $1}')
num_trimmed_pct=$(cat ${outdir_realpath}/trim_stats/example_trim_stats.txt | grep "Input Read" | cut -d ":" -f3 | awk '{print $1, $2}')
num_trimmed=$(echo $num_trimmed_pct | awk '{print $1}')

# get the reference header info
ref_tag=$(cat ${outdir_realpath}/viral_identification/example_KY369892.1_RV_vid.txt | cut -f3)
ref_header=$(cat ${outdir_realpath}/viral_identification/example_KY369892.1_RV_vid.txt | cut -f4)

# get the reference genome (map_all) id, size, coverage
ref_length=$(grep KY369892.1 ${outdir_realpath}/bbmap_covstats_all_ref/example_map_all_bbmap_covstats.tsv | cut -f3)
ref_coverage=$(grep KY369892.1 ${outdir_realpath}/bbmap_covstats_all_ref/example_map_all_bbmap_covstats.tsv | cut -f5)
median_coverage=$(grep KY369892.1 ${outdir_realpath}/bbmap_covstats_all_ref/example_map_all_bbmap_covstats.tsv | cut -f10)

# get the consensus final length
consensus_length=$(awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' ${outdir_realpath}/consensus_final/example_KY369892.1_RV.consensus_final.fa | awk 'FNR==2{print val,$1}')

# get the number and the percentage of trimmed reads mapped to the reference genome (map_ref)
mapped_reads=$(cat ${outdir_realpath}/map_consensus_final_stats/example_KY369892.1_RV_mapf_stats.txt | grep "mapped:" | cut -f3)
# works for single-end and paired-end reads
mapped_reads=$(echo $mapped_reads | awk '{print $1+$2}')

# calculate % reads on target for single-end and paired-end
if [[ false == false ]]
then
percent_mapped_reads=$(echo "$mapped_reads/$num_trimmed*100" | bc -l | awk 'FNR==1{print val,$1}')
else
percent_mapped_reads=$(echo "$mapped_reads/$num_trimmed/2*100" | bc -l | awk 'FNR==1{print val,$1}')
fi

# get the number of Ns in the consensus final 
num_ns=$(grep -v "^>" ${outdir_realpath}/consensus_final/example_KY369892.1_RV.consensus_final.fa | tr -cd N | wc -c | awk 'FNR==1{print val,$1}')

# get the number of Ns in the consensus final 
ref_num_ns=$(grep -v "^>" ${outdir_realpath}/ref_genome/example_KY369892.1_RV.fa | tr -cd N | wc -c | awk 'FNR==1{print val,$1}')

# get the percentage of Ns in the consensus final
percent_n=$(echo "$num_ns/$consensus_length*100" | bc -l | awk 'FNR==1{print val,$1}')

# get the serotype
serotype=$(awk 'FNR==1{print $1}' ${outdir_realpath}/serotype/example_KY369892.1_RV_serotype.txt)

echo "example	$num_untrimmed	$num_trimmed_pct	KY369892.1	$ref_tag	$ref_header	$ref_length	$ref_num_ns	$ref_coverage	$median_coverage	$consensus_length	$mapped_reads	$percent_mapped_reads	$num_ns	$percent_n	$serotype" >> example_KY369892.1_RV_summary.tsv
