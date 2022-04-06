# REVICA

Revica is a reference-based viral consensus genome assembly pipeline written in [Nextflow](https://www.nextflow.io/). For ease of use, Revica can be run in [Docker](https://docs.docker.com/get-docker/). Revica currently supports consensus genome assembly of rhinovirus (RV), human coronavirus (HCOV), human metapneumovirus (HMPV), human respiratory syncytial virus (HRSV), human parainfluenza virus (HPIV), and Measles morbillivirus (MeV). 

## Methods
Revica consists of the following processes:

- Trim reads with Trimmomatic
- Map trimmed reads to a multifasta reference containing complete genomes of supported viruses using BBMap
- Identify best reference(s) for consensus calling based on median coverage (if different viruses are present, the best reference for each virus is identified and consensus genome is assembled for each virus)
- Consensus genome assembly using Samtools and iVar (consensus calling threshold: minimum coverage of 3, minimum base quality of 15, and minimum frequency threshold of 0.6)
- Determine serotypes by BLASTing to a curated BLAST database. 

## Usage
Examples:

	nextflow run greninger-lab/revica -r main --reads input_fastq/fastq.gz_dir_path --outdir output_dir_path


or with Docker


	nextflow run greninger-lab/revica -r main --reads input_fastq/fastq.gz_dir_path --outdir output_dir_path -with-docker greningerlab/revica
	
Valid command line arguments:

	REQUIRED:
	--reads				Input fastq or fastq.gz directory path
	--outdir		        Output directory path 

	OPTIONAL:
	--pe				For paired-end reads (default: single-end)
	--ref				Overwrite reference file
	--m				The median coverage threshold for the initial reference to be considered (default 5)
	--deduplicate			Deduplicated reads using picard before consensus genome assembly
	--help				Displays help message

## Usage notes
- You can use your own reference(s) for consensus genome assembly instead of using our curated database by specifying the `--ref` parameter followed by your fasta file. 
	- it's important to tag the fasta sequences of the same species with the same name or abbreviation in the header section, otherwise the pipeline
	will mistaken references for the the same species but with different accessions as coinfection/multiple infection and generate a consensus
	genome for every reference where the median coverage of the first alignment exceed the specified threshold (default 5).  
	- if you would like to use Revica for a virus not listed, give it a try by supplying a fasta file of reference(s). For segmented viral genomes,
	you could try running individual segments at a time. (Serotype can't be determined for viruses not listed)
- If you are using Docker on Linux, check these [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/) (especially cgroup swap limit capabilities support) for configuring Linux to work better with Docker. 
- By default, Docker has full access to full RAM and CPU resources of the host. However, if you have Docker Desktop installed, go to Settings -> Resources to make sure enough resources (>4 cpus & >4 GB RAM) are allocated to docker containers. 
- More (or less) RAM and CPU sources can be allocated to each process in the pipeline by modifying the `nextflow.config` file.
- Nextflow processes run in separate folders in the `work` directory. This is where you can inspect all the generated files of each process for each input. 
- The `work` directory can take up a lot of space, routinely clean them up if you are done with your analysis. 
- Some Nextflow tips:
	- if you have Revica downloaded locally, you can run the pipeline using the following command
	`nextflow run pipeline_dir_path/main.nf`
	- use `nextflow pull greninger-lab/revica` to download or update to the latest Revica release
	- use `-resume` to resume a run; this requires the `work` directory
	- use `-profile` to select differernt configurations in the `nextflow.config` file
	- use `-r` to specify different github revisions, branches, and releases
	- use `nextflow log <run name> option` to show [information about pipeline execution](https://www.nextflow.io/docs/latest/tracing.html)
	- use `-with-report` to see [resource usage metrics](https://www.nextflow.io/docs/latest/metrics.html)


## Output notes
- The names of the output files are prefixed with `base_refid_refsp` where
	- base: the base name of the input fastq.gz file
	- refid: the accession number of the reference
	- refsp: the abbreviated name of reference species
- Consensus genomes are in the `consensus_final` folder.
- Sorted bam files of trimmed reads to consensus genome are in the `consensus_final_bam_sorted` folder.
- The references used for consensus genome assembly are in the `ref_genome` folder.
- Sorted bam files of trimmed reads to their references are in the `map_ref_bam_sorted` folder.
- Samples with failed assembly are in the `failed_assembly` folder. These samples have less then 5 median coverage to the initial reference. 
- Statistics are in the `run_summary` folder.

## Docker
Pull the docker image (255MB) that has all the dependencies from Docker Hub using the following command

	docker pull greningerlab/revica

Or build the docker image by navigating to the `docker` directory and run the following command

	docker build -t greningerlab/revica .

(You can specify any name after `-t` argument using the format `name:tag`)

To run the docker image interactively, run the following command

	docker run -it --rm greningerlab/revica

[Revica docker image on Docker Hub](https://hub.docker.com/repository/docker/greningerlab/revica) 

## Contact
For bug reports please email aseree@uw.edu or raise an issue on Github.


