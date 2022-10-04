#!/bin/bash

# Map to the multifasta reference
bbmap.sh \
    in=example.trimmed.fastq.gz \
    outm=example_map_all.sam \
    ref=ref.fa \
    threads=8 \
    covstats=example_map_all_bbmap_covstats.tsv \
    local=true interleaved=false maxindel=80 -Xmx8g > example_map_all_stats.txt 2>&1

# sort bbmap_out.txt based on median_fold (column 10) in descending order 
head -1 example_map_all_bbmap_covstats.tsv > example_map_all_bbmap_covstats.temp.txt
awk 'NR>1' < example_map_all_bbmap_covstats.tsv | sort -t $'	' -nrk10 >> example_map_all_bbmap_covstats.temp.txt
mv example_map_all_bbmap_covstats.temp.txt example_map_all_bbmap_covstats.tsv
