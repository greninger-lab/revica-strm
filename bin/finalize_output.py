#!/usr/bin/env python3

import argparse
import os
from collections import defaultdict
import csv
import time
import sys
import subprocess
import re
import glob

SEGMENTED_PATTERN = "H[0-9]{1}N[0-9]{1}"

def collect_segmented_files(sample, files):
    """
    for a given sample, identify all files (usually fasta or bam) 
    that are associated with a segmented reference genome, and return them in a list to be merged 
    into a bam/fasta for the entire organism (all segments)

    just checks the filename for flu
    """
    segment_files = []

    for i, file in enumerate(files):
        if re.search(SEGMENTED_PATTERN, file):
            print(f"segment: {file}")
            segment_files.append(file)
            del files[i]

    return segment_files


def sample_bams_to_fastq(sample, bams, threads):
    """import all mapped reads to a fastq.gz file from all bam files of a given sample"""
    fq_name = f"{sample}.SRA.fastq.gz"

    with open(fq_name, "a") as outf:
        for bam in bams:
            cmd = [
                    "samtools",
                    "fastq",
                    "-F", 
                    "4",
                    "-N",
                    "-@", 
                    threads,
                    bam
                    ]

            print(f"bam to fastq: {cmd}")
            subprocess.check_call(cmd, stdout=outf)

def merge_sample_bams(sample, bams, threads):
    if bams:
        merge_file = f"{sample}_MER_.bam"

        cmd = [
                "samtools",
                "merge",
                "-f",
                merge_file,
                ] + bams

        print(f"merge bams: {cmd}")
        subprocess.run(cmd)
        return merge_file

    return None


def index_bam(bam, threads):
    cmd = [
            "samtools",
            "index",
            bam,
            "-@", 
            threads
            ]

    subprocess.run(cmd)

def match_with_samples(file, samples_files):
    for s in samples_files.keys():
        if file.startswith(s):
            samples_files[s].append(file)


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
def merge_bams(samples, files, threads, merge):
    # samples_bam = defaultdict(list)
    segments = None
    samples_bam = init_sample_dict(samples)
    temp_bams = []
    for file in files:
        if file.endswith(".bam"):
            match_with_samples(file, samples_bam)

    print(f"SAMPLES_MAP: {samples_bam}")
    for sample in samples_bam:
        # don't merge (or remove) files from unsegmented genomes
        if len(samples_bam[sample]) > 1 and merge:
            segments = collect_segmented_files(sample, samples_bam[sample])
            merge_file = merge_sample_bams(sample, segments, threads)
            if merge_file: 
                index_bam(merge_file, threads)
            temp_bams.extend(samples_bam[sample])

        if segments: samples_bam[sample].extend(segments)
        sample_bams_to_fastq(sample, samples_bam[sample], threads)

# combine per-segment final consensus fastas into a multifasta
# for all samples with segmented genomes
def merge_fastas(samples, files, merge, suffix):
    # samples_fa = defaultdict(list)
    samples_fa = init_sample_dict(samples)
    temp_fastas = []
    for file in files:
        if file.endswith(".fa"):
            match_with_samples(file, samples_fa)

    print(f"SAMPLES_MAP: {samples_fa}")
    for sample in samples_fa:
        if len(samples_fa[sample]) > 1 and merge:
            multifasta = os.path.join(
                sample + f"_assembly_MER_{suffix}.fa"
            )

            segments = collect_segmented_files(sample, samples_fa[sample])

            with open(multifasta, "w") as outfile:
                for file in segments:
                    try:
                        with open(file, "r") as readfile:
                            outfile.write(readfile.read())
                        temp_fastas.append(file)
                    except FileNotFoundError:
                        pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("samplesheet")
    parser.add_argument("files", nargs="+")
    parser.add_argument("--suffix", type=str, required=False, default="")
    parser.add_argument("--threads", type=str)
    parser.add_argument("--merge", action="store_true", required=False)

    args = parser.parse_args()

    print("FINALIZING OUTPUT")
    print(args)

    try:
        samples = get_samples(args.samplesheet)
    except FileNotFoundError:
        sys.exit(
            "Samplesheet not found! Check if the samplesheet was moved during processing."
        )

    merge_bams(samples, args.files, args.threads, args.merge)
    merge_fastas(samples, args.files, args.merge, args.suffix)
