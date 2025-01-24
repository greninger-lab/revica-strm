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


def init_sample_dict(samples):
    sample_dict = defaultdict(list)
    for s in samples:
        sample_dict[s] = []

    return sample_dict


def get_samples(samplesheet):
    with open(samplesheet, newline="") as file:
        reader = csv.reader(file)
        next(reader)
        samples = []

        for row in reader:
            samples.append(row[0])  # samples should be first column of samplesheet

    if len(samples) == 0:
        sys.exit(
            "No samples found in first column of samplesheet! Is it appropriately formatted?"
        )

    return samples


# combine the per-segment bam files into a multibam for
# all samples with segmented genomes
def merge_bams(bam_dir, samples):
    # samples_bam = defaultdict(list)
    samples_bam = init_sample_dict(samples)
    temp_bams = []
    for file in os.listdir(bam_dir):
        if file.endswith(".bam") and "_MER_" not in file:
            match_with_samples(bam_dir, file, samples_bam)

    for sample in samples_bam:
        # don't merge (or remove) files from unsegmented genomes
        if len(samples_bam[sample]) > 1:
            merge_file = os.path.join(
                os.path.basename(sample) + "_to_ref_MER_.bam"
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
    # samples_fa = defaultdict(list)
    samples_fa = init_sample_dict(samples)
    temp_fastas = []
    for file in os.listdir(fasta_dir):
        if file.endswith("_final.fa") and "_MER_" not in file:
            match_with_samples(fasta_dir, file, samples_fa)

    for sample in samples_fa:
        if len(samples_fa[sample]) > 1:
            multifasta = os.path.join(
                os.path.basename(sample) + "_assembly_MER_final.fa"
            )

            with open(multifasta, "w") as outfile:
                for file in samples_fa[sample]:
                    try:
                        with open(file, "r") as readfile:
                            outfile.write(readfile.read())
                        temp_fastas.append(file)
                        print(file)
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
        samples = get_samples(args.samplesheet)
    except FileNotFoundError:
        sys.exit(
            "Samplesheet not found! Check if the samplesheet was moved during processing."
        )

    merge_bams(BAM_DIR, samples)
    merge_fastas(FASTA_DIR, samples)
