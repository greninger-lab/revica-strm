#!/bin/bash -ue
echo "example" > sample_name.list

bcftools mpileup \
    --fasta-ref example_KY369892.1_RV.fa \
    --ignore-overlaps \
    --count-orphans \
    --no-BAQ \
    --max-depth 0 \
    --min-BQ 20 \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    example_KY369892.1_RV_map_ref.sorted.bam \
    | bcftools call \
        --output-type v \
        --ploidy 1 \
        --keep-alts \
        --keep-masked-ref \
        --multiallelic-caller \
        --variants-only \
    | bcftools reheader \
        --samples sample_name.list \
    | bcftools view \
        --output-file example.vcf.gz \
        --output-type z \
        --include 'INFO/DP>=10'

tabix -p vcf -f example.vcf.gz

bcftools stats example.vcf.gz > example.bcftools_stats.txt
