# REVICA-STRM

STRM = Streamlined

---

Revica is a reference-based viral consensus genome assembly pipeline for some of the most common respiratory viruses. Revica currently supports genome assembly of:
- Enterovirus (EV)
- Seasonal human coronavirus (HCOV)
- Human metapneumovirus (HMPV)
- Human respiratory syncytial virus (HRSV)
- Human parainfluenza virus (HPIV)
- Measles morbillivirus (MeV)
- Influenza A virus (IAV)
- Influenza B virus (IBV)
- Human adenovirus (HAdV)

## How it works

REVICA-STRM creates assembly genomes from raw FASTQ files in 2 fundamental steps:

1. **Create an initial, rough consensus sequence.** A given query sequence is aligned to all entries in a given database, with the best-mapping reference used to create an early consensus sequence with the general predicted features.

2. **Generate the final consensus for the query**. The query is realigned to the consensus from the previous step, and this alignment is used to generate the final assembly. This second pass serves to confirm unique features of the query that may have been suppressed during the previous assembly.
---
The important outputs are found in the output directory (specified by `--output`), in the `final_files` subdirectory:
- the final consensus sequence (filename contains `assembly.fa`)
- the BAM file(s) generated from the first alignment: this is useful for gauging how much of your RAW FASTQ data could be mapped to a known reference sequence. Located in the `align_to_selected_ref` directory.
- a MULTIQC report file named `<run name>_multiqc.html`

## Databases
This tool includes two example reference databases usable for assembly:

#### general respiratory viruses
- `assets/ref.fa`: a general purpose database containing sequences for the aforementioned supported viruses

#### influenza only
- `assets/flu.fasta`: the database used in the [Andersen Lab's avian influenza project](https://github.com/andersen-lab/avian-influenza). This database has been supplemented with recent influenza A H5, H1N1 and H3N1 strains, as well as Influenza B, across all available species collected from 2023-2024.

## Workflow
![Workflow](revica_workflow_diagram.png)

## Usage

REVICA is built to be run on the cloud via NextFlow and Docker. Cloning this repo is only necessary if you want the databases, scripts, or test data.

- Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- Install [`Docker`](https://docs.docker.com/engine/installation/)
- **it's recommended to use Docker signed out**; access to certain containers is sporadically blocked if signed in. This issue is being actively investigated.
- Ensure the Docker client is running before starting the pipeline

### Using the test data included in this repo:


1. Clone the repository to get the example data and database  
    `git clone git@github.com:epiliper/nf-rev.git`

    `cd nf-rev`
2. Run the pipeline with the example data
    ```bash
    nextflow run greninger-lab/REVICA-STRM -r main -latest --input example_samplesheet.csv --output example_output -profile docker --db assets/flu.fasta 
    ```

After the run has finished, the final output files can be found in `<work_folder, default=run>/final_files`. 

If not using example data, replace the FASTQ files, sample sheet, and database with whatever files you want to use (**see below**).

### Using other FASTQ files and databases:

Cloning this repo is not necessary unless you need the example data. 

1. download FASTQ files for needed samples/SRA projects. The [SRA toolkit's](https://github.com/ncbi/sra-tools) `fasterq_dump` utility can be used for downloading FASTQ files from SRA projects.

2. in the directory with the downloaded FASTQs, use this repo's included script `bin/fastq_dir_to_samplesheet.py` to create a REVICA sample sheet.

    Example command:
    ```bash
    python3 fastq_dir_to_samplesheet.py <dir with fastqs> -r1 _1.fastq.gz -r2 _2.fastq.gz sras_to_run.csv     
    ```
    
    where `-r1` and `-r2` specify the suffixes of input read1 and read2 files to look for, respectively. For single-end data, just use `-r1`.

3. run REVICA and point it to your sample sheet:

    ```bash
    nextflow run greninger-lab/REVICA-STRM -r main -latest --input sras_to_run.csv --output example_output -profile docker --db assets/flu.fasta
    ```

>[!Note]  
>This repo includes a python script `bin/pull_sra.py` to download FASTQ files for SRA project numbers specified in a CSV spreadsheet, and create an associated REVICA spreadsheet.    
>
>To use it, download SRA toolkit, add it to $PATH, and populate a csv file with SRAs in the format according to `assets/example_sras.csv`
>
>then run `python3 pull_sra.py <sra csv> <name for REVICA samplesheet>`
## Options
|Option|Explanation|
|------|-----------|
| `--input` | samplesheet in csv format with fastq information |
| `--output` | output directory (default: revica_output) |
| `--db` | (multi)fasta file to overwrite the bundled viral database |
| `--run_name` | name for the summary tsv file (default: 'run') |
| `--skip_fastp` | skip adapters and reads trimming using fastp (default: false) |
| `--run_kraken2` | run Kraken2 for classifying reads (default: false) |
| `--kraken2_db` | Kraken2 database for reads classification, needs to be specified when using `--run_kraken2` |
| `--kraken2_variants_host_filter` | use reads that didn't map to the kraken2 database for downstream consensus calling |
| `--save_kraken2_unclassified_reads` | save reads that didn't map to the specified kraken2 database |
| `--save_kraken2_classified_reads` | save reads that map to the specified kraken2 database |
| `--trim_len` | minimum read length to keep (default:50) |
| `--save_trimmed_reads` | save trimmed fastq |
| `--save_temp_files` | save temporary files |
| `--sample` | downsample fastq to a certain fraction or number of reads |
| `--ref_min_median_cov` | minimum median coverage on a reference for consensus assembly (default: 3) |
| `--ref_min_genome_cov` | minimum reference coverage percentage for consensus assembly (default: 60%) |
| `--ivar_consensus_t` | minimum frequency threshold to call consensus (default: 0.6) |
| `--ivar_consensus_q` | minimum quality score threshold to call consensus (default: 20) |
| `--ivar_consensus_m` | minimum depth to call consensus (default: 5) |

## Usage notes
- Samplesheet example: `assets/samplesheet.csv`
- You can create a samplesheet using the bundled python script: `python bin/fastq_dir_samplesheet.py fastq_dir samplesheet_name.csv`
- Memory and CPU usage for pipeline processes can be adjusted in `conf/base.config`
- Process arguments can be adjusted in `conf/modules.config`
- You can use your own reference(s) for consensus genome assembly by specifying the `--db` parameter followed by your fasta file. 
	- reference header format: `>reference_accession reference_tag reference_header_info`
	- it's important to tag the fasta sequences for the same species or gene segments with the same name or abbreviation in the header section, otherwise the pipeline
	will generate a consensus genome for every reference where the median coverage of the first alignment exceed the specified threshold (default 3).  
	- Revica works with segmented viral genomes, just keep the different gene segments separated and tag them in the reference fasta file
- If you are using Docker on Linux, check out these [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/) (especially cgroup swap limit capabilities support) for configuring Linux to work better with Docker. 
- By default, Docker has full access to full RAM and CPU resources of the host, but if you are using MacOS, go to Settings -> Resources in Docker Desktop to make sure enough resources are allocated to docker containers. 

## Contact
For bug reports, please raise an issue.
