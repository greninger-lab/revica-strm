import subprocess
import argparse
import os
import gzip
import shutil
import pandas as pd

# downloads the fastqs from the SRAs in the csv
# compresses them to fastq.gz for compatibility with revica

parser = argparse.ArgumentParser()
parser.add_argument("sra")
parser.add_argument("name")

args = parser.parse_args()


def pull_sra(sra):
    print(f"Downloading {sra} fastq(s)...", end="\r")
    subprocess.run(["fasterq-dump", "--format", "fastq", sra], capture_output=True)
    print("compressing files to .fastq.gz...", end="\r")
    for file in os.listdir("."):
        if file.startswith(sra) and file.endswith(".fastq"):
            outfile = file.split(".")[0] + ".fastq.gz"
            with open(file, "rb") as f_in:
                with gzip.open(outfile, "wb", compresslevel=1) as f_out:
                    shutil.copyfileobj(f_in, f_out)

            os.remove(file)


# pull the sra fastqs
df = pd.read_csv(args.sra)
num_sras = len(df["SRA"])

for i, sra in enumerate(df["SRA"], 1):
    print(f"processing entry {i} of {num_sras}:")
    print("------------------------------------")
    pull_sra(sra)
    print("\033[2A\r", end="\r")

# create a REVICA samplesheet from the pulled, zipped fastqs
subprocess.run(
    [
        "python3",
        "fastq_dir_to_samplesheet.py",
        ".",
        "-r1",
        "_1.fastq.gz",
        "-r2",
        "_2.fastq.gz",
        # "sras_to_run.csv",
        args.name,
    ]
)
