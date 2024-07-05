#!/usr/bin/env python

import argparse
import os
import subprocess

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")

    args = parser.parse_args()

    # outdir = os.path.relpath(os.getcwd(), os.path.join(args.outdir))
    samples_bam = {}
    samples_fa = {}

    bam_dir = os.path.join(args.outdir, "bbmap_align", "consensus")
    fasta_dir = os.path.join(args.outdir, "ivar_consensus")

    # merge the bam files for each segment into a multibam
    for file in os.listdir(bam_dir):
        if file.endswith(".bam") and "merged" not in file:
            prefix = file.split("_")[0]
            segment = os.path.join(bam_dir, file)

            if prefix not in samples_bam:
                samples_bam[prefix] = [segment]
            else:
                samples_bam[prefix].append(segment)

    temp_bams = []
    temp_fastas = []

    for sample in samples_bam:
        # don't merge (or remove) files from unsegmented genomes
        if len(samples_bam[sample]) > 1:
            merge_file = os.path.join(bam_dir, f"{sample}_merged.bam")
            subprocess.run(
                [
                    "samtools",
                    "merge",
                    "-f",
                    merge_file,
                    *samples_bam[sample],
                ],
                shell=True,
            )
            # pysam.merge("-@ 4", "-f", merge_file, *samples_bam[sample])
            temp_bams.append(merge_file)

            # [os.remove(segment_bam) for segment_bam in samples_bam[sample]]

    # remove indexes from segment bams
    [
        os.remove(os.path.join(bam_dir, index))
        for index in os.listdir(bam_dir)
        if index.endswith(".bai")
    ]

    for file in os.listdir(fasta_dir):
        if file.endswith(".fa"):
            prefix = file.split("_")[0]
            segment = os.path.join(fasta_dir, file)

            if prefix not in samples_fa:
                samples_fa[prefix] = [segment]

            else:
                samples_fa[prefix].append(segment)

    for sample in samples_fa:
        if len(samples_fa[sample]) > 1:
            multifasta = os.path.join(fasta_dir, sample + "_merged.fa")

            with open(multifasta, "w") as outfile:
                for file in samples_fa[sample]:
                    with open(file, "r") as readfile:
                        outfile.write(readfile.read())

                    temp_fastas.append(file)

    [os.remove(segment_fa) for segment_fa in temp_fastas]
