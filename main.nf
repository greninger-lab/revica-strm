#!/usr/bin/env nextflow

/*
========================================================================================
                    		   REVICA v1.0
========================================================================================
Github Repo:
https://github.com/greninger-lab/revica

Author:
Jaydee Sereewit <aseree@uw.edu>
Alex L Greninger <agrening@uw.edu>
UW Medicine | Virology
Department of Laboratory Medicine and Pathology
University of Washington
Updated: March 30, 2022
Revision: v1.0	
LICENSE: GNU
----------------------------------------------------------------------------------------
*/
// Pipeline version
version = '1.0'
// Setup pipeline help msg
def help() {
    log.info"""
    _______________________________________________________________________________
    				REVICA :  Version ${version}
    _______________________________________________________________________________
    
	Pipeline Usage:
	Natively
	nextflow run greninger-lab/revica -r main --reads input_fastq.gz_path --outdir output_dir_path
	or with Docker (easy but slower than running natively):
	nextflow run greninger-lab/revica -r main --reads input_fastq.gz_path --outdir output_dir_path -with-docker greningerlab/revica


    Valid CLI Arguments:
    REQUIRED:
	--reads			Input fastq.gz path
	--outdir		Output directory path
    OPTIONAL:
	--pe			Specifies that the input fastq files are paired-end reads
	--deduplicate		Get rid of duplicated reads before consensus genome assembly
	--ref_rv		Overwrite set multifasta rhinovirus reference file
	--ref_hcov		Overwrite set multifasta human coronavirus reference file
	--ref_hmpv		Overwrite set multifasta human metapneumovirus reference file
	--ref_hrsv		Overwrite set multifasta human respiratory syncytial virus reference file
	--ref_hpiv		Overwrite set multifasta human parainfluenza virus reference file
	--help			Displays help message in terminal
    """
}

// Setup command-line parameters
params.help = false
params.reads = false
params.outdir = false
params.pe = false
params.deduplicate = false
params.ref_rv = false
params.ref_hcov = false
params.ref_hmpv = false
params.ref_hrsv = false
params.ref_hpiv = false

// Check Nextflow version for enabling DSL2
nextflow_dsl2_v = '20.07.1'
if ( nextflow.version.matches(">= $nextflow_dsl2_v") ) {
    nextflow.enable.dsl=2
} else {
    nextflow.preview.dsl=2
}

// Show help msg if --help parameter is set or if reads AND outdir are not specified
if (params.help == true || (params.reads == false && params.outdir == false)){
    help()
    exit 0
}

// if INPUT not set
if (params.reads == false) {
    println( "Please provide input fastq.gz path with --reads") 
    exit(1)
}

// if OUTDIR not set
if (params.outdir == false) {
    println( "Please provide output directory path with --outdir") 
    exit(1)
}

// Setup MULTIFASTA Reference file paths for OPTIONAL-override of set file paths.
// Rhinovirus
if(params.ref_rv != false) {
    Reference_rv = file(params.ref_rv)
} else {
Reference_rv=file("${baseDir}/ref/hrv_ref_rv.fa")
}
// Human Coronavirus
if(params.ref_hcov != false) {
    Reference_hcov = file(params.ref_hcov)
} else {
Reference_hcov=file("${baseDir}/ref/hrv_ref_hcov.fa")
}
// Human metapneumovirus 
if(params.ref_hmpv != false) {
    Reference_hmpv = file(params.ref_hmpv)
} else {
Reference_hmpv=file("${baseDir}/ref/hrv_ref_hmpv.fa")
}
// Human respiratory syncytial virus 
if(params.ref_hrsv != false) {
    Reference_hrsv = file(params.ref_hrsv)
} else {
Reference_hrsv=file("${baseDir}/ref/hrv_ref_hrsv.fa")
}
// Human parainfluenza virus 
if(params.ref_hpiv != false) {
    Reference_hpiv = file(params.ref_hpiv)
} else {
Reference_hpiv=file("${baseDir}/ref/hrv_ref_hpiv.fa")
}
// Adapters file path
ADAPTERS = file("${baseDir}/adapters/adapters.fa")

//Trimmomatic Settings
params.SETTING = "2:30:10:1:true"
params.LEADING = "3"
params.TRAILING = "3"
params.SWINDOW = "4:20"
params.MINLEN = "35"

// bbmap maxindel for different viruses
params.rv_maxindel = 9
params.hcov_maxindel = 18
params.hrsv_maxindel = 60
params.hmpv_maxindel = 40
params.default_maxindel = 80

// Setup Blast database for serotyping
// All BLAST db files for respiratory viruses recognized by this pipeline including:
// rhinovirus
// human coronavirus 229E, OC43, NL63, HKU1
// human respiratory syncytial virus A, B
// human metapneumovirus A1, A2, B1, B2 
BLAST_DB_ALL_1 = file("${baseDir}/blast_db/all_ref.fasta")
BLAST_DB_ALL_2 = file("${baseDir}/blast_db/all_ref.fasta.ndb")
BLAST_DB_ALL_3 = file("${baseDir}/blast_db/all_ref.fasta.nhr")
BLAST_DB_ALL_4 = file("${baseDir}/blast_db/all_ref.fasta.nin")
BLAST_DB_ALL_5 = file("${baseDir}/blast_db/all_ref.fasta.not")
BLAST_DB_ALL_6 = file("${baseDir}/blast_db/all_ref.fasta.nsq")
BLAST_DB_ALL_7 = file("${baseDir}/blast_db/all_ref.fasta.ntf")    
BLAST_DB_ALL_8 = file("${baseDir}/blast_db/all_ref.fasta.nto")  

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
include { Serotyping } from './modules.nf'
include { Summary_Generation } from './modules.nf'
include { Final_Processing } from './modules.nf'

// Create channel for input reads: single-end or paired-end
if(params.pe == false) {
    // Looks for gzipped files, assumes all separate samples
    input_read_ch = Channel
        .fromPath("${params.reads}/*.gz")
        .ifEmpty { error "Cannot find any FASTQ in ${params.reads} ending with .gz" }
        .map { it -> file(it)}
} else {
    // Check for R1s and R2s in input directory
    input_read_ch = Channel
        .fromFilePairs("${params.reads}/*_R{1,2}*.gz")
        .ifEmpty { error "Cannot find any FASTQ pairs in ${params.reads} ending with .gz" }
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
            ADAPTERS,
            params.MINLEN,
            params.SETTING, 
            params.LEADING,
            params.TRAILING,
            params.SWINDOW
        )

        Aligning_SE (
            Trimming_SE.out[0],
            Reference_rv,
            Reference_hcov, 
	    Reference_hmpv,
	    Reference_hrsv,
	    Reference_hpiv
	)

	Viral_Identification (
	    Aligning_SE.out[0],
	    Aligning_SE.out[1]
	)

	Consensus_Generation_Prep_SE (
	    Viral_Identification.out[0].flatten(),
	    Reference_rv,
	    Reference_hcov,
	    Reference_hmpv,
	    Reference_hrsv,
	    Reference_hpiv
	)

        Consensus_Generation_SE (Consensus_Generation_Prep_SE.out[0])

	Serotyping (
	    Consensus_Generation_SE.out[0],
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

        Final_Processing (Summary_Generation.out[0].collect())

    } else { 
	
	Trimming_PE (
	    input_read_ch, 
            ADAPTERS,
            params.MINLEN,
            params.SETTING, 
            params.LEADING,
            params.TRAILING,
            params.SWINDOW
	)

        Aligning_PE (
            Trimming_PE.out[0],
            Reference_rv,
            Reference_hcov, 
	    Reference_hmpv,
	    Reference_hrsv,
	    Reference_hpiv
	)

	Viral_Identification (
	    Aligning_PE.out[0],
	    Aligning_PE.out[1]
	)

	Consensus_Generation_Prep_PE (
	    Viral_Identification.out[0].flatten(),
	    Reference_rv,
	    Reference_hcov,
	    Reference_hmpv,
	    Reference_hrsv,
	    Reference_hpiv
	)

        Consensus_Generation_PE (Consensus_Generation_Prep_PE.out[0])

	Serotyping (
	    Consensus_Generation_PE.out[0],
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

        Final_Processing (Summary_Generation.out[0].collect())

    }
}
