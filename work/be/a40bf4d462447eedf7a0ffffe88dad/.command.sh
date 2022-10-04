#!/bin/bash

base=$(head -1 *_vid.txt | cut -f1)
ref_id=$(head -1 *_vid.txt | cut -f2)
ref_tag=$(head -1 *_vid.txt | cut -f3)
samtools faidx ref.fa $ref_id > ${base}'_'${ref_id}'_'${ref_tag}'.fa'

# capture process environment
set +u
echo base=${base[@]} > .command.env
echo ref_id=${ref_id[@]} >> .command.env
echo ref_tag=${ref_tag[@]} >> .command.env
