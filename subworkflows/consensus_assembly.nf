//                                                                                 
// 3 iterations of consensus assembly using BBMAP for alignment and iVar consensus for consensus calling
//                                                              

include { BBMAP_ALIGN as BBMAP_ALIGN_QUERY                              } from '../modules/bbmap_align'
include { IVAR_CONSENSUS_BBMAP_ALIGN as IVAR_CONSENSUS_BBMAP_ALIGN_INITIAL_ASSEMBLY        } from './ivar_consensus_bbmap_align'                          
include {IVAR_CONSENSUS as BUILD_FINAL_CONSENSUS } from '../modules/ivar_consensus'
//include { IVAR_CONSENSUS_BBMAP_ALIGN as IVAR_CONSENSUS_BBMAP_ALIGN_2        } from './ivar_consensus_bbmap_align'                          
// include { IVAR_CONSENSUS_BBMAP_ALIGN as IVAR_CONSENSUS_BBMAP_ALIGN_FINAL_ASSEMBLY    } from './ivar_consensus_bbmap_align'                          
                                                                                   
workflow CONSENSUS_ASSEMBLY {                                                 
    take:                                                                          
    ch_reads    // channel: [ val(meta), path(reads) ]                  
    ch_ref      // channel: [ val(meta), val(ref_info), path(ref) ]                             
                                                                                
    main:                                                                          

    BBMAP_ALIGN_QUERY (
        ch_reads,
        ch_ref
    )

    IVAR_CONSENSUS_BBMAP_ALIGN_INITIAL_ASSEMBLY (
        BBMAP_ALIGN_QUERY.out.bam,
        BBMAP_ALIGN_QUERY.out.ref,
        BBMAP_ALIGN_QUERY.out.reads
    )

 //   IVAR_CONSENSUS_BBMAP_ALIGN_FINAL_ASSEMBLY (
 //       IVAR_CONSENSUS_BBMAP_ALIGN_INITIAL_ASSEMBLY.out.bam,
 //       IVAR_CONSENSUS_BBMAP_ALIGN_INITIAL_ASSEMBLY.out.consensus,
 //       IVAR_CONSENSUS_BBMAP_ALIGN_INITIAL_ASSEMBLY.out.reads
 //   )

    BUILD_FINAL_CONSENSUS (
        IVAR_CONSENSUS_BBMAP_ALIGN_INITIAL_ASSEMBLY.out.bam,
        IVAR_CONSENSUS_BBMAP_ALIGN_INITIAL_ASSEMBLY.out.consensus,
    )

    emit:
    // consensus   = IVAR_CONSENSUS_BBMAP_ALIGN_INITIAL_ASSEMBLY.out.consensus // channel: [ val(meta), val(ref_info), path(consensus) ]
    consensus   = BUILD_FINAL_CONSENSUS.out.consensus
    bam         = IVAR_CONSENSUS_BBMAP_ALIGN_INITIAL_ASSEMBLY.out.bam       // channel: [ val(meta), val(ref_info), path(bam), path(bai) ]
    reads       = IVAR_CONSENSUS_BBMAP_ALIGN_INITIAL_ASSEMBLY.out.reads     // channel: [ val(meta), path(reads) ]
} 
