#!/usr/bin/env python
import argparse
import os
import pysam
from collections import defaultdict
import csv
import sys


def match_with_samples(dir, file, samples_files):
    for s in samples_files.keys():
        if file.startswith(s):
            samples_files[s].append(os.path.join(dir, file))


def get_samples(csv_reader_samplesheet):
    next(csv_reader_samplesheet)
    samples = []

    for row in csv_reader_samplesheet:
        samples.append(row[0])  # samples should be first column of samplesheet

    if len(samples) == 0:
        sys.exit(
            "No samples found in first column of samplesheet! Is it appropriately formatted?"
        )


# combine the per-segment bam files into a multibam for
# all samples with segmented genomes
def merge_bams(bam_dir, samples):
    samples_bam = defaultdict(list)
    temp_bams = []
    for file in os.listdir(bam_dir):
        if file.endswith(".bam") and "_MER_" not in file:
            match_with_samples(bam_dir, file, samples_bam)

    for sample in samples_bam:
        # don't merge (or remove) files from unsegmented genomes
        if len(samples_bam[sample]) > 1:
            merge_file = os.path.join(
                bam_dir, os.path.basename(sample) + "_to_ref_MER_.bam"
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
def merge_fastas(fasta_dir, samples):
    samples_fa = defaultdict(list)
    temp_fastas = []
    for file in os.listdir(fasta_dir):
        if file.endswith("_final.fa") and "_MER_" not in file:
            match_with_samples(fasta_dir, file, samples_fa)

    for sample in samples_fa:
        if len(samples_fa[sample]) > 1:
            multifasta = os.path.join(
                fasta_dir, os.path.basename(sample) + "_assembly_MER_final.fa"
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
    parser.add_argument("samplesheet")

    args = parser.parse_args()

    BAM_DIR = os.path.join(args.outdir, "final_files", "align_to_selected_ref")
    FASTA_DIR = os.path.join(args.outdir, "ivar_consensus")

    try:
        reader = csv.reader(args.samplesheet)
    except FileNotFoundError:
        sys.exit(
            "Samplesheet not found! Check if the samplesheet was moved during processing."
        )

    samples = get_samples(reader)

    merge_bams(BAM_DIR, samples)
    merge_fastas(FASTA_DIR, samples)
