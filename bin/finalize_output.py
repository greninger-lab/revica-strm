#!/usr/bin/env python3

import argparse
import os
from collections import defaultdict
import csv
import time
import sys
import subprocess

def bam_to_fastq(files, threads):
    for file in files:
        fq_name = os.path.basename(file).split(".bam")[0] + ".SRA.fastq.gz"

        cmd = [
                "samtools",
                "fastq",
                "-F", 
                "4",
                file,
                "-o",
                fq_name,
                "-N",
                "-@", 
                threads
                ]

        print(f"bam to fastq: {cmd}")
        subprocess.run(cmd, stderr=subprocess.DEVNULL)

def merge_sample_bams(sample, bams, threads):
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

def index_bam(bam, threads):
    cmd = [
            "samtools",
            "index",
            bam,
            "-@", 
            threads
            ]

    subprocess.run(cmd)

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
def merge_bams(bam_dir, samples, threads, merge):
    # samples_bam = defaultdict(list)
    samples_bam = init_sample_dict(samples)
    temp_bams = []
    for file in os.listdir(bam_dir):
        if file.endswith(".bam") and "_MER_" not in file:
            match_with_samples(bam_dir, file, samples_bam)

    for sample in samples_bam:
        # don't merge (or remove) files from unsegmented genomes
        if len(samples_bam[sample]) > 1 and merge:
            merge_file = merge_sample_bams(sample, samples_bam[sample], threads)
            index_bam(merge_file, threads)
            bam_to_fastq([merge_file], threads)
            temp_bams.extend(samples_bam[sample])
        else: 
            bam_to_fastq(samples_bam[sample], threads)

# combine per-segment final consensus fastas into a multifasta
# for all samples with segmented genomes
def merge_fastas(fasta_dir, samples, merge):
    # samples_fa = defaultdict(list)
    samples_fa = init_sample_dict(samples)
    temp_fastas = []
    for file in os.listdir(fasta_dir):
        if file.endswith("_final.fa") and "_MER_" not in file:
            match_with_samples(fasta_dir, file, samples_fa)

    for sample in samples_fa:
        if len(samples_fa[sample]) > 1 and merge:
            multifasta = os.path.join(
                os.path.basename(sample) + "_assembly_MER_final.fa"
            )

            with open(multifasta, "w") as outfile:
                for file in samples_fa[sample]:
                    try:
                        with open(file, "r") as readfile:
                            outfile.write(readfile.read())
                        temp_fastas.append(file)
                    except FileNotFoundError:
                        pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")
    parser.add_argument("samplesheet")
    parser.add_argument("threads", type=str)
    parser.add_argument("--merge", action="store_true", required=False)

    args = parser.parse_args()

    print("FINALIZING OUTPUT")
    print(args)

    # this directory should have reads aligned to the initial assemblies, and 
    # the initial assembly fastas themselves
    BAM_DIR = os.path.join(args.outdir, "final_files", "align_to_consensus")
    FASTA_DIR = os.path.join(args.outdir, "final_files", "final_assemblies")

    try:
        samples = get_samples(args.samplesheet)
    except FileNotFoundError:
        sys.exit(
            "Samplesheet not found! Check if the samplesheet was moved during processing."
        )

    # merge all BAM files aligned to initial consensus, and 
    # merge the initial assemblies themselves for use as BAM reference
    merge_bams(BAM_DIR, samples, args.threads, args.merge)
    merge_fastas(BAM_DIR, samples, args.merge)

    merge_fastas(FASTA_DIR, samples, args.merge)
