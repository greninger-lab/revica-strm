#!/usr/bin/env python
import argparse
import os
import pysam


def get_file_prefix(file):
    if os.path.abspath(file):
        return os.path.basename(file).split("_")[0]

    else:
        return file.split("_")[0]


# combine the per-segment bam files into a multibam for
# all samples with segmented genomes
def merge_bams(bam_dir):
    samples_bam = {}
    for file in os.listdir(bam_dir):
        if file.endswith(".bam") and "merged" not in file:
            # prefix = file.split("_")[0]
            prefix = get_file_prefix(file)
            segment = os.path.join(bam_dir, file)

            if prefix not in samples_bam:
                samples_bam[prefix] = [segment]
            else:
                samples_bam[prefix].append(segment)

    for sample in samples_bam:
        # don't merge (or remove) files from unsegmented genomes
        if len(samples_bam[sample]) > 1:
            merge_file = os.path.join(bam_dir, os.path.basename(sample) + "_merged.bam")
            # subprocess.run(
            #     [
            #         "samtools",
            #         "merge",
            #         "-f",
            #         "-o",
            #         merge_file,
            #     ]
            #     + samples_bam[sample],
            #     shell=True,
            #     check=True,
            # )
            pysam.merge("-f", merge_file, *samples_bam[sample])

    for sample in samples_bam:
        for temp in samples_bam[sample]:
            try:
                os.remove(temp)
            except FileNotFoundError:
                pass


# combine per-segment final consensus fastas into a multifasta
# for all samples with segmented genomes
def merge_fastas(fasta_dir):
    samples_fa = {}
    temp_fastas = []
    for file in os.listdir(fasta_dir):
        if file.endswith(".fa") and "merged" not in file:
            prefix = file.split("_")[0]
            segment = os.path.join(fasta_dir, file)

            if prefix not in samples_fa:
                samples_fa[prefix] = [segment]

            else:
                samples_fa[prefix].append(segment)

    for sample in samples_fa:
        if len(samples_fa[sample]) > 1:
            multifasta = os.path.join(
                fasta_dir, os.path.basename(sample) + "_merged.fa"
            )

            with open(multifasta, "w") as outfile:
                for file in samples_fa[sample]:
                    try:
                        with open(file, "r") as readfile:
                            outfile.write(readfile.read())
                    except FileNotFoundError:
                        pass
                    temp_fastas.append(file)

    for temp in temp_fastas:
        try:
            os.remove(temp)
        except FileNotFoundError:
            pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")

    args = parser.parse_args()

    bam_dir = os.path.join(args.outdir, "bbmap_align", "consensus")
    fasta_dir = os.path.join(args.outdir, "ivar_consensus")

    merge_bams(bam_dir)
    merge_fastas(fasta_dir)
