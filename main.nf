#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-BIRD
    Version: 1.0.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/Kincekara/nfbird
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CBIRD  } from './workflows/cbird'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {    
    
    main:
    if (params.read1 && params.read2 && params.samplename){
        // Create a channel directly from the provided parameters
        ch_samplesheet = Channel.of([params.samplename, file(params.read1), file(params.read2)])        
    } else if (params.samplesheet){
        // Parse the samplesheet CSV file to create the channel
        Channel.fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> 
            [row.sample, [file(row.fastq_1), file(row.fastq_2)] ] 
        }
        .set { ch_samplesheet }               
    } else {
        error "No input file provided. Please provide a samplesheet file." 
    }
    //
    // WORKFLOW: Run main workflow
    //
    CBIRD (ch_samplesheet)
}

