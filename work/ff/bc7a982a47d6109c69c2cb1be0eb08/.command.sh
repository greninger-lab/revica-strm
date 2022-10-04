#!/bin/bash

# get total input reads for mapping
reads_used=$(grep "Reads Used:" example_map_all_stats.txt | cut -f2)

# get mapped reads
mapped_reads=$(grep "Mapped reads:" example_map_all_stats.txt | cut -f2)

# identify initial reference
python3 /Users/greningerlab/AS/revica/bin/select_reference.py -bbmap_covstats example_map_all_bbmap_covstats.tsv -b example -reads_num ${reads_used} -mapped_reads ${mapped_reads} -m 3 -p 0
