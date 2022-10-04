#!/usr/bin/env nextflow

/*
========================================================================================
                    		   REVICA v1.2.2
========================================================================================
Github Repo:
https://github.com/greninger-lab/revica

Author:
Jaydee Sereewit <aseree@uw.edu>
Alex L Greninger <agrening@uw.edu>
UW Medicine | Virology
Department of Laboratory Medicine and Pathology
University of Washington
Updated: Septemper 29, 2022
LICENSE: GNU
----------------------------------------------------------------------------------------
*/
// Pipeline version
version = '1.2.2'
// Setup pipeline help msg
def help() {
    log.info"""
    _______________________________________________________________________________
    				REVICA :  Version ${version}
    _______________________________________________________________________________
    
	Pipeline Usage:
	Natively
	nextflow run greninger-lab/revica -r main -latest --reads input_fastq/fastq.gz_dir_path --outdir output_dir_path
	or with Docker:
	nextflow run greninger-lab/revica -r main -latest --reads input_fastq/fastq.gz_dir_path --outdir output_dir_path -with-docker greningerlab/revica


    Valid CLI Arguments:
    REQUIRED:
	--reads			Input fastq or fastq.gz directory path
	--outdir		Output directory path
    OPTIONAL:
	--pe			For paired-end reads (default: single-end)
	--ref			Overwrite reference file
	--dedup         	Get rid of duplicated reads before consensus genome assembly
	--sample		Subsample reads to a fraction or a number
	--m			The median coverage threshold for the initial reference to be considered (default: 3)
	--p                     The minimum covered percent by the reads for the initial refernce to be considered (default: 0)
        --q                     Minimum base quality score threshold for iVar consensus to count base. (Default: 20)
        --t                     Minimum frequency threshold(0 - 1) for iVar call consensus. (Default: 0.6)
        --d                     Minimum depth for iVar to call consensus. (Default: 5)
	--help			Displays help message in terminal
    """
}

// Setup command-line parameters
params.reads = null
params.outdir = null
params.pe = false
params.ref = false
params.dedup = false
params.sample = false
params.m = 3
params.p = 0
params.q = 20
params.t = 0.6
params.d = 5
params.help = false

//Trimmomatic Settings
params.SETTING = "2:30:10:1:true"
params.LEADING = "3"
params.TRAILING = "3"
params.SWINDOW = "4:20"
params.MINLEN = 35

// Check Nextflow version for enabling DSL2
nextflow_dsl2_v = '20.07.1'
if ( nextflow.version.matches(">= $nextflow_dsl2_v") ) {
    nextflow.enable.dsl=2
} else {
    nextflow.preview.dsl=2
}

// Set color to red on default background
def fg = 31
def bg = 49
def color = "${(char)27}[$fg;$bg"+"m"

// Show help msg if --help parameter is set or if reads AND outdir are not specified
if (params.help == true || (params.reads == false && params.outdir == false)){
    help()
    exit 0
}

// if INPUT not set
if (params.reads == null) {
    println(color+"Please provide input fastq.gz path with --reads") 
    exit(1)
}

// if OUTDIR not set
if (params.outdir == null) {
    println(color+"Please provide output directory path with --outdir") 
    exit(1)
}

// Setup MULTIFASTA Reference file path for OPTIONAL-override of set file path.
if (params.ref != false) {
    ref = file(params.ref)
} else {
    ref = file("${baseDir}/ref/ref.fa")
}

// Trimmomatic adapters file path
ADAPTERS = file("${baseDir}/adapters/adapters.fa")

// Setup Blast database for rhinovirus serotyping
BLAST_DB_ALL_1 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta")
BLAST_DB_ALL_2 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta.ndb")
BLAST_DB_ALL_3 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta.nhr")
BLAST_DB_ALL_4 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta.nin")
BLAST_DB_ALL_5 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta.not")
BLAST_DB_ALL_6 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta.nsq")
BLAST_DB_ALL_7 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta.ntf")    
BLAST_DB_ALL_8 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta.nto")  

// Setup pipeline header
def header() {
    return """
    """.stripIndent()
}
// log files header
// log.info header()
log.info "________________________________________________________________________________"
log.info "                             REVICA :  v${version}"
log.info "________________________________________________________________________________"
def summary = [:]
summary['Nextflow run name'] = workflow.runName
summary['Revica directory'] = workflow.projectDir
if(workflow.revision) summary['Pipeline release'] = workflow.revision
summary['Configuration profile'] = workflow.profile
summary['Launch directory'] = workflow.launchDir
summary['Input directory'] = params.reads
summary['Sequence type'] = params.pe ? 'Paired-End' : 'Single-End'
summary['Output directory'] = params.outdir
summary['Work directory'] = workflow.workDir
summary['Trimmomatic adapters'] = ADAPTERS
summary["Trimmomatic min read length"] = params.MINLEN
summary["Trimmomatic setting"] = params.SETTING
summary["Trimmomatic sliding Window"] = params.SWINDOW
summary["Trimmomatic leading"] = params.LEADING
summary["Trimmomatic trailing"] = params.TRAILING
log.info summary.collect { k,v -> "${k.padRight(30)}: $v" }.join("\n")
log.info "________________________________________________________________________________"

// Import processes
include { Trimming_SE } from './modules.nf'
include { Trimming_PE } from './modules.nf'
include { Aligning_SE } from './modules.nf'
include { Aligning_PE } from './modules.nf'
include { Viral_Identification } from './modules.nf'
include { Consensus_Generation_Prep_SE } from './modules.nf'
include { Consensus_Generation_Prep_PE } from './modules.nf'
include { Consensus_Generation_SE } from './modules.nf'
include { Consensus_Generation_PE } from './modules.nf'
include { VCF_Generation } from './modules.nf'
include { Serotyping } from './modules.nf'
include { Summary_Generation } from './modules.nf'
include { Final_Processing } from './modules.nf'

// Create channel for input reads: single-end or paired-end
if(params.pe == false) {
    // Looks for gzipped files, assumes all separate samples
    input_read_ch = Channel
        .fromPath("${params.reads}/*{.fastq.gz,.fastq}")
        .ifEmpty { error "Cannot find any FASTQ in ${params.reads} ending with fastq or fastq.gz" }
        .map { it -> file(it)}
} else {
    // Check for R1s and R2s in input directory
    input_read_ch = Channel
        .fromFilePairs("${params.reads}/*_{R1,R2,1,2}*{.fastq.gz,.fastq}")
        .ifEmpty { error "Cannot find any FASTQ pairs in ${params.reads} ending with fastq or fastq.gz" }
        .map { it -> [it[0], it[1][0], it[1][1]]}
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    if(params.pe == false) {

	Trimming_SE (
            input_read_ch, 
            ADAPTERS
        )

        Aligning_SE (
            Trimming_SE.out[0],
	    ref
	)

	Viral_Identification (
	    Aligning_SE.out[0],
	    ref
	)

	Consensus_Generation_Prep_SE (
	    Viral_Identification.out[0].flatten(),
	    ref
	)

        Consensus_Generation_SE (
            Consensus_Generation_Prep_SE.out[0]
        )

        VCF_Generation (
            Consensus_Generation_SE.out[0]
        )

	Serotyping (
	    Consensus_Generation_SE.out[1],
	    BLAST_DB_ALL_1,
	    BLAST_DB_ALL_2,
	    BLAST_DB_ALL_3,
	    BLAST_DB_ALL_4,
	    BLAST_DB_ALL_5,
	    BLAST_DB_ALL_6,
	    BLAST_DB_ALL_7,
	    BLAST_DB_ALL_8
	)

        Summary_Generation (
	    Serotyping.out[0]
	)

        Final_Processing (
            Summary_Generation.out[0].collect()
        )

    } else { 
	
	Trimming_PE (
	    input_read_ch, 
            ADAPTERS
	)

        Aligning_PE (
            Trimming_PE.out[0],
	    ref
	)

	Viral_Identification (
	    Aligning_PE.out[0],
	    ref
	)

	Consensus_Generation_Prep_PE (
	    Viral_Identification.out[0].flatten(),
	    ref
	)

        Consensus_Generation_PE (
            Consensus_Generation_Prep_PE.out[0]
        )

        VCF_Generation (
            Consensus_Generation_PE.out[0]
        )

	Serotyping (
	    Consensus_Generation_PE.out[1],
	    BLAST_DB_ALL_1,
	    BLAST_DB_ALL_2,
	    BLAST_DB_ALL_3,
	    BLAST_DB_ALL_4,
	    BLAST_DB_ALL_5,
	    BLAST_DB_ALL_6,
	    BLAST_DB_ALL_7,
	    BLAST_DB_ALL_8
	)

        Summary_Generation (
	    Serotyping.out[0]
	)

        Final_Processing (
            Summary_Generation.out[0].collect()
        )

    }
}
