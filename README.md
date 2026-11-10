# revica-strm

strm = streamlined

This is a fork of the [original REVICA](https://github.com/greninger-lab/revica), meant for faster assembly of viral genomes from short-read sequencing data. It uses the BWA MEM aligner and employs additional assembly quality checks to provide higher-quality consensus genomes, while striving to make output files and reporting more concise and informative.

**Update 2025 July 16**: BWA MEM is now the default aligner used, due to [this issue](https://github.com/bwa-mem2/bwa-mem2/issues/262) with BWA MEM2. To use BWA MEM2, pass the `--use_mem2` option.

---

revica-strm is a reference-based viral consensus genome assembly pipeline for some of the most common respiratory viruses. revica-strm default references datasets currently supports genome assembly of:

- Enterovirus (EV, RV)
- Parechovirus (HPeV)
- Seasonal human coronavirus (229E, HKU1, NL63, OC43)
- SARS-CoV-2
- Human metapneumovirus (HMPV)
- Human respiratory syncytial virus (HRSV)
- Human parainfluenza virus (HPIV1-4)
- Measles morbillivirus (MeV)
- Mumps Virus (MuV)
- Rubella virus (RuV)
- Influenza A virus (H1N1, H3N2, H5N1, H5N5, H7N9)
- Influenza B virus (Victoria, Yamagata)
- Influenza C virus
- Human adenovirus (HAdV A-G)
- Human Bocavirus
- Anellovirus
- Herpesvirus (HSV1, HSV2, HHV-6, HHV-7, CMV, EBV, VZV, KSHV)
- Polyomavirus (WUPyV, KIPyV)
- Nipah virus (NiV)

## How it works

revica-strm creates assembly genomes from raw FASTQ files in 2 fundamental steps:

1. **Create an initial, rough consensus sequence.** A given query sequence is aligned to all entries in a given database, with the best-mapping reference used to create an early consensus sequence with the general predicted features.

2. **Generate the final consensus for the query**. The query is realigned to the consensus from the previous step, and this alignment is used to generate the final assembly. This second pass serves to confirm unique features of the query that may have been suppressed during the previous assembly.
---
The important outputs are found in the output directory (specified by `--output`), in the `final_files` subdirectory:
- the final consensus sequence (filename contains `assembly.fa`)
- the BAM file(s) generated from the first alignment: this is useful for gauging how much of your RAW FASTQ data could be mapped to a known reference sequence. Located in the `align_to_selected_ref` directory.
- a MULTIQC report file named `<run name>_multiqc.html`

## Databases
This repository includes two example reference databases usable for assembly:

#### general respiratory viruses
- `assets/ref_20251119_reheadered.fa`: Our most up-to-date database with a balance between sequence redundance and species diversity. It is the default.

#### influenza only
- `assets/flu.fasta`: the database used in the [Andersen Lab's avian influenza project](https://github.com/andersen-lab/avian-influenza). This database has been supplemented with recent influenza A H5, H1N1 and H3N1 strains, as well as Influenza B, across all available species collected from 2023-2024.

## Workflow
![Workflow](revica_workflow_diagram.png)

## Usage

`revica-strm <DIRECTORY WITH FASTQS> -profile docker [OPTIONS]`

Run `revica-strm -h` for option details.

## Installation

1. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) if you haven't already
2. Install [Docker](https://docs.docker.com/desktop/) if you haven't already
3. download the latest `revica-strm` run script and make it executable by running this command:

```bash
wget https://raw.githubusercontent.com/greninger-lab/revica-strm/refs/heads/main/revica-strm
chmod +x revica-strm
revica-strm
```

**For convenience, it's recommended to move `revica-strm` to somewhere in your `$PATH`, so you can run it from other directories**. Otherwise, you'll have to run it with something like `eval path/to/revica-strm`. For instructions on adding to `$PATH`, see [here](https://superuser.com/questions/488173/how-can-i-edit-the-path-on-linux).


## Instructions

1. Ensure the docker desktop client is updated and running

2. Arrange all input FASTQs (can be plain .fastq or compressed .fastq.gz) in their own directory.

##### 🚨 MANDATORY: all FASTQ files must have unique sample names before the first underscore ('\_') character.

RIGHT: ```sample1_R1.fastq.gz sample1_R2.fastq.gz```  
WRONG: ```sample_1_R1.fastq.gz sample_2_R1.fastq.gz```  

- In the wrong case, sample_1_R1 would get wrongly paired with sample_2_R1.
- *It's recommended to have the read mate info (e.g. R1, 1) immediately follow the first underscore, after the unique sample name.*
- if sample names do not have underscores, this logic applies to the period ('.') character instead.

3. once you're sure your FASTQs are correctly named, run:

```
revica-strm <fastq_dir> -profile docker --output <out_dir>
```

You can test this command by downloading the example dataset included in this repository, `fastq_example`, and running:
```
revica-strm fastq_example -profile docker
```

For a description of all options, run `revica-strm -h`.

Once the pipeline is done, consensus genomes can be found in `$out_dir/final_files/final_assemblies/`. For reports of any genomes/samples that failed assembly QC thresholds, see files in `$out_dir/fail`.

## Notes on storage use
Nextflow tracks its internal progress via the creation of a work directory (default = `./work`). This can grow to be quite large; it is advised that if you are operating with limited disk space to batch your datasets into smaller runs and delete the work directory in-between. Note that the work directory enables Nextflow to run from the last error/crash-point with the `-resume` flag, so the work directory should only be deleted if you are certain that you are done with your dataset and/or have no more space.

## 2026-Nov-10: checking for redundant consensus sequences
When using a reference DB with entries that are highly similar but have distinct ref tags (e.g. H5N1_PB1 and H5N5_PB1), REVICA-STRM can occasionally assemble multiple consensus sequences for a sample that really has just one organism. REVICA-STRM now has an optional flag (`--check_dup`) that will cluster sequences by their dedup tag (see [Reference database](#Reference-database)), perform pairwise alignments on each cluster, and report consensus sequences per-sample that have high pairwise-identity (default: >= 99.9%) in `$out_dir/dedup`. Within each cluster, consensus sequences are ranked by coverage, then depth; each is grouped under the first higher-ranked consensus it matches, and that higher-ranked consensus is reported as the distinct consensus. 

This is a purely post-processing step with no impact on upstream steps: all output consensus sequences are preserved regardless of similarity, pending review using the dedup report.

## 🤓 For developers

### Running the pipeline with advanced nextflow options

The `revica-strm` script is just a wrapper over the basic `nextflow run` command; it creates a fastq samplesheet from given input directory, and passes other arguments to `nextflow run` itself. You can pass any of the typical `nextflow run` arguments to the `revica-strm` script. This includes the `-c` command to specify advanced options for, as an example, running the pipeline on AWS Batch or other cloud computing environments.

### Removing host (human) reads
Inputs to revica-strm can optionally be filtered with Kraken2 and a user-supplied Kraken2 database. This database should be comprised of host/contaminant genomes desired to be removed from downstream analysis.

To use this, run revica-strm with the `--run-kraken2` and `--kraken2_variants_host_filter` commands, and point the `--kraken2_db` argument to your kraken2 database.

>[!NOTE]
>To create a database we recommend for removal of human reads, see [these instructions](making_kraken2_human_db.md).

### Reference database
The reference database is comprised of multiple representatives of a variety of respiratory virus species such as enterovirus, seasonal coronavirus, SARS-CoV2, parainfluenza, measles, influenza, and more. Inspect `assets/ref_20251119_reheadered.fa` if curious. If you intend to use your own database, ensure the fasta headers are structured as follows:   

```ACCESSION<SPACE>DEDUP_TAG<SPACE>REF_TAG<SPACE>SAMPLE_HEADER<SPACE>|<SPACE>NCBI_TITLE```   

where REF_TAG should be unique to a species-specific segment/genome, and DEDUP_TAG groups entries that should be checked against each other for redundant consensus sequences (see [checking for redundant consensus sequences](#2026-Nov-10-checking-for-redundant-consensus-sequences)). We use `<species_taxid>[_<segment>]` (e.g. all Flu A PB1 entries share one DEDUP_TAG, regardless of subtype). SAMPLE_HEADER is what is reported in the `ref_header` column of the run summary. The `| NCBI_TITLE` part is optional and is ignored by the pipeline. Take these entries for Flu A segments PB1 and NS1, and an enterovirus genome, for example:

```
>KX023793.1 2955291_PB1 FluA_H1N1_PB1 A/DC/2/2016 | Influenza A virus (A/DC/2/2016(H1N1)) segment 2 polymerase PB1 (PB1) gene, complete cds; and nonfunctional PB1-F2 protein (PB1-F2) gene, complete sequence
>KX023799.1 2955291_NS FluA_H1N1_NS A/DC/2/2016 | Influenza A virus (A/DC/2/2016(H1N1)) segment 8 nuclear export protein (NEP) and nonstructural protein 1 (NS1) genes, complete cds
>PQ895337.1 3428506 EV_D68 SH2024-25870 | Enterovirus D68 isolate SH2024-25870, complete genome
```

#### Converting an existing database
`bin/reheader_ref_fasta.py` converts a database in the old `ACCESSION<SPACE>REF_TAG<SPACE>SAMPLE_HEADER` format to the format above. It looks up each accession's species taxid (and the record title) from NCBI, builds the DEDUP_TAG (the segment is parsed from the old REF_TAG), and keeps the old tag and header as REF_TAG and SAMPLE_HEADER. It also writes a lookup table of all the old and new values:

```
bin/reheader_ref_fasta.py -i ref.fa -o ref.reheadered.fa -t ref.lookup.tsv [-k $NCBI_API_KEY]
```

You can specify your own reference database with `--refs $REF_DB`.


## Contact
For bug reports, please raise an issue or contact me via the links in my [profile](https://github.com/epiliper/):
