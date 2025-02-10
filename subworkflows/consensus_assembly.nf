//                                                                                 
// 2 iterations of consensus assembly using BBMAP for alignment and iVar consensus for consensus calling
//                                                              

include { BWA_MEM2_ALIGN as BWA_MEM2_ALIGN_QUERY                                     } from '../modules/bwa_mem2_align' 
include { IVAR_CONSENSUS_BWA_ALIGN as IVAR_CONSENSUS_BWA_MEM2_ALIGN_INITIAL_ASSEMBLY } from './ivar_consensus_bwa_align'
include { IVAR_CONSENSUS as BUILD_FINAL_CONSENSUS                                    } from '../modules/ivar_consensus'
                                                                                   
workflow CONSENSUS_ASSEMBLY {                                                 
    take:                                                                          
    ch_reads    // channel: [ val(meta), path(reads) ]                  
    ch_ref      // channel: [ val(meta), val(ref_info), path(ref) ]                             
                                                                                
    main:                                                                          

    BWA_MEM2_ALIGN_QUERY(
        ch_reads,
        ch_ref
    )

    IVAR_CONSENSUS_BWA_MEM2_ALIGN_INITIAL_ASSEMBLY(
        BWA_MEM2_ALIGN_QUERY.out.bam,
        BWA_MEM2_ALIGN_QUERY.out.ref,
        BWA_MEM2_ALIGN_QUERY.out.reads
    )


    BUILD_FINAL_CONSENSUS (
        IVAR_CONSENSUS_BWA_MEM2_ALIGN_INITIAL_ASSEMBLY.out.bam,
        IVAR_CONSENSUS_BWA_MEM2_ALIGN_INITIAL_ASSEMBLY.out.consensus
    )

    emit:
    final_consensus     = BUILD_FINAL_CONSENSUS.out.consensus
    initial_consensus   = IVAR_CONSENSUS_BWA_MEM2_ALIGN_INITIAL_ASSEMBLY.out.consensus
    bam                 = IVAR_CONSENSUS_BWA_MEM2_ALIGN_INITIAL_ASSEMBLY.out.bam       // channel: [ val(meta), val(ref_info), path(bam), path(bai) ]
    reads               = IVAR_CONSENSUS_BWA_MEM2_ALIGN_INITIAL_ASSEMBLY.out.reads     // channel: [ val(meta), path(reads) ]
} 
