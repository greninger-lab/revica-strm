//                                                                                 
// Consensus Assembly using BBMAP for alignment and iVar consensus for consensus calling  
//                                                              
                                                                                   
include { IVAR_CONSENSUS    } from '../modules/ivar_consensus'                          
include { BBMAP_ALIGN       } from '../modules/bbmap_align' 
                                                                                   
workflow IVAR_CONSENSUS_BBMAP_ALIGN {                                                 
    take:                                                                          
    ch_bam                // channel: [ val(meta), val(ref_info), path(bam), path(bai) ]                             
    ch_ref                // channel: [ val(meta), val(ref_info), path(ref) ]                             
    ch_reads              // channel: [ val(meta), path(reads) ]                  
                                                                                   
    main:                                                                          

    IVAR_CONSENSUS (
        ch_bam,
        ch_ref
    )

    ch_reads
        .join(IVAR_CONSENSUS.out.consensus)
        .multiMap { meta, reads, ref_info, consensus -> 
            reads:      [ meta, reads ]
            consensus:  [ meta, ref_info, consensus]
        }
        .set { ch_bbmap_align_input }

    BBMAP_ALIGN (
        ch_bbmap_align_input.reads,
        ch_bbmap_align_input.consensus
    )
    
    emit:
    bam         = BBMAP_ALIGN.out.bam       // channel: [ val(meta), val(ref_info), path(bam), path(bai) ]
    consensus   = BBMAP_ALIGN.out.ref       // channel: [ val(meta), val(ref_info), path(consensus) ]
    reads       = BBMAP_ALIGN.out.reads     // channel: [ val(meta), path(reads) ]
}
