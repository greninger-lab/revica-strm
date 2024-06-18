//                                                                                 
// 3 iterations of consensus assembly using BBMAP for alignment and iVar consensus for consensus calling
//                                                              

include { BBMAP_ALIGN as BBMAP_ALIGN_REFERENCE                              } from '../modules/bbmap_align'
include { IVAR_CONSENSUS_BBMAP_ALIGN as IVAR_CONSENSUS_BBMAP_ALIGN_1        } from './ivar_consensus_bbmap_align'                          
include { IVAR_CONSENSUS_BBMAP_ALIGN as IVAR_CONSENSUS_BBMAP_ALIGN_2        } from './ivar_consensus_bbmap_align'                          
include { IVAR_CONSENSUS_BBMAP_ALIGN as IVAR_CONSENSUS_BBMAP_ALIGN_FINAL    } from './ivar_consensus_bbmap_align'                          
                                                                                   
workflow CONSENSUS_ASSEMBLY {                                                 
    take:                                                                          
    ch_reads    // channel: [ val(meta), path(reads) ]                  
    ch_ref      // channel: [ val(meta), val(ref_info), path(ref) ]                             
                                                                                
    main:                                                                          

    BBMAP_ALIGN_REFERENCE (
        ch_reads,
        ch_ref
    )

    IVAR_CONSENSUS_BBMAP_ALIGN_1 (
        BBMAP_ALIGN_REFERENCE.out.bam,
        BBMAP_ALIGN_REFERENCE.out.ref,
        BBMAP_ALIGN_REFERENCE.out.reads
    )

    IVAR_CONSENSUS_BBMAP_ALIGN_2 (
        IVAR_CONSENSUS_BBMAP_ALIGN_1.out.bam,
        IVAR_CONSENSUS_BBMAP_ALIGN_1.out.consensus,
        IVAR_CONSENSUS_BBMAP_ALIGN_1.out.reads
    )

//   IVAR_CONSENSUS_BBMAP_ALIGN_FINAL (
//       IVAR_CONSENSUS_BBMAP_ALIGN_2.out.bam,
//       IVAR_CONSENSUS_BBMAP_ALIGN_2.out.consensus,
//       IVAR_CONSENSUS_BBMAP_ALIGN_2.out.reads
//   )

    emit:
    consensus   = IVAR_CONSENSUS_BBMAP_ALIGN_1.out.consensus // channel: [ val(meta), val(ref_info), path(consensus) ]
    bam         = IVAR_CONSENSUS_BBMAP_ALIGN_1.out.bam       // channel: [ val(meta), val(ref_info), path(bam), path(bai) ]
    reads       = IVAR_CONSENSUS_BBMAP_ALIGN_1.out.reads     // channel: [ val(meta), path(reads) ]
} 
