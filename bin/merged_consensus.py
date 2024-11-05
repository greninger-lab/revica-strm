#!/usr/bin/env python
import argparse
import os
import pysam
from collections import defaultdict


def get_file_prefix(file):
    return os.path.basename(file).split("_")[0]


# combine the per-segment bam files into a multibam for
# all samples with segmented genomes
def merge_bams(bam_dir):
    samples_bam = defaultdict(list)
    temp_bams = []
    for file in os.listdir(bam_dir):
        if file.endswith(".bam") and "merged" not in file:
            prefix = get_file_prefix(file)
            segment = os.path.join(bam_dir, file)

            samples_bam[prefix].append(segment)

    for sample in samples_bam:
        # don't merge (or remove) files from unsegmented genomes
        if len(samples_bam[sample]) > 1:
            merge_file = os.path.join(
                bam_dir, os.path.basename(sample) + "_to_ref_merged.bam"
            )
            pysam.merge("-f", merge_file, *samples_bam[sample])
            temp_bams.extend(samples_bam[sample])

    for temp in temp_bams:
        try:
            os.remove(temp)
        except FileNotFoundError:
            pass


# combine per-segment final consensus fastas into a multifasta
# for all samples with segmented genomes
def merge_fastas(fasta_dir):
    samples_fa = defaultdict(list)
    temp_fastas = []
    for file in os.listdir(fasta_dir):
        if file.endswith("_final.fa") and "merged" not in file:
            prefix = get_file_prefix(file)
            segment = os.path.join(fasta_dir, file)

            samples_fa[prefix].append(segment)

    for sample in samples_fa:
        if len(samples_fa[sample]) > 1:
            multifasta = os.path.join(
                fasta_dir, os.path.basename(sample) + "_assembly_merged_final.fa"
            )

            with open(multifasta, "w") as outfile:
                for file in samples_fa[sample]:
                    try:
                        with open(file, "r") as readfile:
                            outfile.write(readfile.read())
                        temp_fastas.append(file)
                    except FileNotFoundError:
                        pass

    for temp in temp_fastas:
        try:
            os.remove(temp)
        except FileNotFoundError:
            pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")

    args = parser.parse_args()

    BAM_DIR = os.path.join(args.outdir, "final_files", "align_to_selected_ref")
    FASTA_DIR = os.path.join(args.outdir, "ivar_consensus")

    merge_bams(BAM_DIR)
    merge_fastas(FASTA_DIR)
