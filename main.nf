#!/usr/bin/env nextflow
/*
   ========================================================================================
   REVICA
   ========================================================================================
   Github Repo:
https://github.com/asereewit/revica

Author:
Jaydee Sereewit <aseree@uw.edu>
Alex L Greninger <agrening@uw.edu>
UW Medicine | Virology
Department of Laboratory Medicine and Pathology
University of Washington
LICENSE: GNU
----------------------------------------------------------------------------------------
 */

// if INPUT not set
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// if db not set

if (!params.db) { exit 1, "Reference database not specified!"}
//
// SUBWORKFLOWS
//
include { INPUT_CHECK               } from './subworkflows/input_check'
include { REFERENCE_PREP            } from './subworkflows/reference_prep'
include { CONSENSUS_ASSEMBLY        } from './subworkflows/consensus_assembly'
include { FINALIZE_OUTPUT           } from './modules/finalize_output'

//
// MODULES
//
include { SEQTK_SAMPLE              } from './modules/seqtk_sample'
include { SUMMARY                   } from './modules/summary'
include { KRAKEN2                   } from './modules/kraken2'
include { FASTQ_TRIM_FASTP_MULTIQC  } from './subworkflows/fastq_trim_fastp_multiqc.nf'
include { DELETE_TEMP_FILES         } from './modules/delete_temp_files'

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
log.info "                                                                                        "
log.info " ██████╗ ███████╗██╗   ██╗██╗ ██████╗ █████╗       ███████╗████████╗██████╗ ███╗   ███╗ "
log.info " ██╔══██╗██╔════╝██║   ██║██║██╔════╝██╔══██╗      ██╔════╝╚══██╔══╝██╔══██╗████╗ ████║ "
log.info " ██████╔╝█████╗  ██║   ██║██║██║     ███████║█████╗███████╗   ██║   ██████╔╝██╔████╔██║ "
log.info " ██╔══██╗██╔══╝  ╚██╗ ██╔╝██║██║     ██╔══██║╚════╝╚════██║   ██║   ██╔══██╗██║╚██╔╝██║ "
log.info " ██║  ██║███████╗ ╚████╔╝ ██║╚██████╗██║  ██║      ███████║   ██║   ██║  ██║██║ ╚═╝ ██║ "
log.info " ╚═╝  ╚═╝╚══════╝  ╚═══╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝ "
log.info "                                                                                        "



workflow {

        INPUT_CHECK (
                ch_input
                )

        FASTQ_TRIM_FASTP_MULTIQC (
                INPUT_CHECK.out.reads,
                params.adapter_fasta,
                params.save_trimmed_fail,
                params.save_merged,
                params.skip_fastp,
                params.skip_fastqc
                )

        ch_sample_input = FASTQ_TRIM_FASTP_MULTIQC.out.reads

        if (params.run_kraken2) {
            KRAKEN2 (
                    FASTQ_TRIM_FASTP_MULTIQC.out.reads,
                    file(params.kraken2_db),
                    params.kraken2_variants_host_filter || params.kraken2_assembly_host_filter,
                    params.kraken2_variants_host_filter || params.kraken2_assembly_host_filter
                    )

                if (params.kraken2_variants_host_filter) {
                    ch_sample_input = KRAKEN2.out.unclassified_reads_fastq
                }
        }

    if (params.sample) {
        SEQTK_SAMPLE (
                ch_sample_input,
                params.sample
                )
            ch_ref_prep_input = SEQTK_SAMPLE.out.reads
    } else {
        ch_ref_prep_input = ch_sample_input
    } 

    if (!params.skip_consensus) { 
        REFERENCE_PREP (
                ch_ref_prep_input,
                file(params.db)
                ) 

            CONSENSUS_ASSEMBLY (
                    REFERENCE_PREP.out.reads,
                    REFERENCE_PREP.out.ref,
                    )

            FASTQ_TRIM_FASTP_MULTIQC.out.trim_log
            .combine(CONSENSUS_ASSEMBLY.out.final_consensus
                    .join(CONSENSUS_ASSEMBLY.out.bam, by: [0,1]), by: 0)
            .map { meta, trim_log, ref_info, consensus, bam, bai -> [ meta, ref_info, trim_log, consensus, bam, bai ] }
        .set { ch_summary_in }

        SUMMARY (
                ch_summary_in
                )

            SUMMARY.out.summary
            .collectFile(storeDir: "${params.output}", name:"${params.run_name}_summary.tsv", keepHeader: true, sort: true)

            SUMMARY.out.ready_to_concat
            .collect()
            .map { it -> true }
        .set { all_summaries_done }

        FINALIZE_OUTPUT(
            CONSENSUS_ASSEMBLY.out.final_consensus.map { meta, ref_info, ref -> [ ref ] }.collect(),
            CONSENSUS_ASSEMBLY.out.initial_consensus.map { meta, ref_info, ref -> [ ref ] }.collect(),
            CONSENSUS_ASSEMBLY.out.bam.map { meta, ref_info, bam, bai -> [ bam ] }.collect(),
            // .map { meta, final_ref_info, final_ref, init_ref_info, init_ref, bam , bai -> [ meta, final_ref, init_ref, bam ] }
            ch_input,
            // file("${params.output}").toAbsolutePath().toString(),
            // file("${params.input}").toAbsolutePath().toString(),
            all_summaries_done,
            params.concat_flu
        )

        if (!params.save_temp_files) {
            DELETE_TEMP_FILES(
                FINALIZE_OUTPUT.out.done,
                file("${params.output}").toAbsolutePath()
            )
        }
    }
}
