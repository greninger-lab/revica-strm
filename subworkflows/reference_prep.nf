//                                                                                 
// Take the db alignment covstats from bbmap and output 
// the selected reference(s) based on selection criteria
// specified in nextflow.config. Copies of reads are
// created to match the number of selected reference(s)
// for downstream alignment and consensus assembly processes.
//

include { BBMAP_ALIGN_DB        } from '../modules/bbmap_align_db'         
include { SELECT_REFERENCE      } from '../modules/select_reference'                     
include { MAKE_REFERENCE_FASTA  } from '../modules/make_reference_fasta'


workflow REFERENCE_PREP {                                                 
    take:                                                                          
    ch_reads    // channel: [ val(meta), path(reads) ]                  
    db          // path: db
                                                                                   
    main:
                                                                                   
    BBMAP_ALIGN_DB (                                                       
        ch_reads,                                                          
        db                                                                 
    )                                                                              
                                                                                   
    SELECT_REFERENCE (                                                             
        BBMAP_ALIGN_DB.out.covstats                                        
    )                                                                              

    // Unpack and reformat the list so each item emitted by the channel is
    // [ [ meta.id, meta.single_end ], [ ref_inf.acc, ref_info.tag, ref_info.header ] ]              
    SELECT_REFERENCE.out.refs_tsv
        .map { meta, refs_tsv -> add_ref_info_to_meta(meta, refs_tsv) }                                                
        .flatten().collate( 2,false )
        .set { ch_new_meta }

    MAKE_REFERENCE_FASTA (                                                                     
        ch_new_meta,                                             
        db                                                                 
    )                                                                              

    // Create copies of reads to match the number of references selected        
    ch_reads
        .cross(MAKE_REFERENCE_FASTA.out.ref)
        .multiMap { it -> 
            reads: it[0]
            ref:   it[1]
        }
        .set { ch_output }

    emit:
    reads   = ch_output.reads   // channel: [ val(meta), path(reads) ]
    ref     = ch_output.ref     // channel: [ val(meta), val(ref_info), path(ref_fasta) ]
}

//
// Parse the output refs_tsv file from SELECT_REFERENCE and return a LIST of
// [ [ meta.id, meta.single_end ], [ ref_inf.acc, ref_info.tag, ref_info.header ] ]
// for the selected reference(s)          
def add_ref_info_to_meta(meta, refs_tsv) {
    def new_meta_list = []

    refs_tsv.text.eachLine { line ->
        def new_meta = []
        def fields = line.split('\t')
        new_meta = [ meta, [ acc: fields[0], tag: fields[1], header: fields[2].toString() ] ]
        new_meta_list.add(new_meta)
    }
    
    return new_meta_list
}
