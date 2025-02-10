#!/usr/bin/env python3

import argparse
import sys
import csv

PANDEPTH_COLS = ["chr", "length", "covered_site", "depth", "coverage", 	"meandepth"]

def parse_pandepth(file, ref, out):

    ref_map = {}

    with open(ref) as refs:
        for line in refs.readlines():
            if line.startswith(">"):
                ref_name, ref_acc = line[1:].strip("\n"), line[1:].strip("\n").split(" ")[0]
                ref_map[ref_acc] = ref_name

    try: 
        with open(file, "r") as inf:
            with open(out, "w") as outf:
                csv_reader = csv.DictReader(inf.readlines()[1:-1], fieldnames=PANDEPTH_COLS, delimiter = "\t")
                csv_writer = csv.DictWriter(outf, fieldnames=PANDEPTH_COLS, delimiter = "\t")

                csv_writer.writeheader()

                for row in csv_reader:
                    # replace just the accesion with the entire name
                    # a miss should never happen unless something REALLY GOES WRONG, 
                    # but will error handle just in case
                    try:
                        row["chr"] = ref_map[row["chr"]]
                    except KeyError:
                        sys.exit("Detected discordance between reference names in pandepth output file and reference! Check work dir output!")

                    csv_writer.writerow(row)

    except FileNotFoundError: 
        sys.exit("Failed to open pandepth file!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help = "unzipped pandepth output file")
    parser.add_argument("ref", help = "reference fasta to use for replacing header")
    parser.add_argument("out", help = "name of output")
    args = parser.parse_args()

    parse_pandepth(args.input, args.ref, args.out)
