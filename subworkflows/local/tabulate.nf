workflow TABULATE {

    take:
    ch_fastp
    ch_assembly_prep
    ch_profile
    ch_quast
    ch_checkm2
    ch_predict_taxon
    ch_mlst
    ch_amrfinder
    ch_plasmidfinder    
    ch_tblastn
    ch_generate_report
    ch_qc_check

    main:
    ch_table = ch_fastp
        .join(ch_assembly_prep, remainder: true)
        .join(ch_profile, remainder: true)
        .join(ch_quast, remainder: true)
        .join(ch_checkm2, remainder: true)
        .join(ch_predict_taxon, remainder: true)
        .join(ch_mlst, remainder: true)
        .join(ch_amrfinder, remainder: true)
        .join(ch_plasmidfinder, remainder: true)
        .join(ch_tblastn, remainder: true)
        .join(ch_generate_report, remainder: true)
        .join(ch_qc_check, remainder: true)

    header = Channel.of(
        "sample_id\tr1_reads\tr2_reads\ttotal_reads\ttotal_reads_trim\tr1_q30_raw\tR1_q30_trim\tr2_q30_raw\tr2_q30_trim\t" +
        "phix_ratio\t" +
        "bracken_taxon\tbracken_taxon_ratio\t" +
        "genome_length\tnumber_of_contigs\tn50_value\tl90_value\tgc_content\t" +
        "completeness\tcontamination\t" +
        "predicted_organism\tpercent_identity\t" +
        "mlst\tpubmlst_scheme\t" +
        "amr_genes\tamr_stress_genes\tamr_virulence_genes\tamr_subclass\t" +
        "plasmidfinder_plasmids\t" +
        "blast_genes\t" +
        "est_genome_length_ratio\test_sequencing_depth\test_sequencing_depth_trim\t" +
        "qc_eval"
        )

    // Combine the header with the data
    header
        .concat(
            ch_table.map { it -> 
                it.collect { value -> value ?: 'NA' } // Replace null with 'NA'
                .join('\t') // Convert each tuple to a tab-separated string
            }
        )
        .collectFile(
            storeDir: "${params.outdir}", 
            name: "results.tsv",               
            sort: false,                       
            newLine: true                       
        )


}

