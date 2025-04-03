process QC_CHECK {
    label 'process_single'
    container 'ubuntu:jammy-20240911.1'
    
    input:
    tuple val(samplename),
        val (genome_ratio),
        val (coverage),
        val (coverage_trim),
        val (completeness),
        val (contamination),
        val (number_of_scaffolds),        
        val (total_reads_trim),
        val (r1_q30_trim),
        val (r2_q30_trim)

    output:
    tuple val(samplename), env('qc_eval'), emit: stats

    script:
    """
    # coverage
    if awk "BEGIN {exit !(${coverage} < 40)}"; then
        echo -n "FAIL:Coverage<40X, " >> QC_EVAL
    elif awk "BEGIN {exit !(${coverage_trim} < 30)}"; then
        echo -n "WARN:Trimmed Coverage<30X, " >> QC_EVAL
    fi
    # contamination
    if awk "BEGIN {exit !(${contamination} > 2.5)}"; then
        echo -n"FAIL:contamination>2.5%, " >> QC_EVAL
    elif awk "BEGIN {exit !(${contamination} > 1.2)}"; then
        echo -n "WARN:contamination>1.2%, " >> QC_EVAL
    fi
    # completeness
    if awk "BEGIN {exit !(${completeness} < 95)}"; then
        echo -n "FAIL:genome_completeness<95%, " >> QC_EVAL
    elif awk "BEGIN {exit !(${completeness} < 97.9)}"; then
        echo -n "WARN:genome_completeness<97.9%, " >> QC_EVAL
    fi
    # q30
    if awk "BEGIN {exit !(${r1_q30_trim} < 90)}"; then
        echo -n "FAIL:R1 Q30<90%, " >> QC_EVAL
    fi
    if awk "BEGIN {exit !(${r2_q30_trim} < 70)}"; then
        echo -n "FAIL:R2 Q30<70%, " >> QC_EVAL
    fi
    # reads
    if awk "BEGIN {exit !(${total_reads_trim} < 1000000)}"; then
        echo -n "WARN:trimmed reads<1M, " >> QC_EVAL
    fi
    # contigs
    if awk "BEGIN {exit !(${number_of_scaffolds} > 200)}"; then
        echo -n "FAIL:Contigs>200, " | >> QC_EVAL
    fi
    # genome ratio
    if awk "BEGIN {exit !(${genome_ratio} > 1.25)}"; then 
        echo -n "WARN:genome ratio>1.25, " >> QC_EVAL
    elif awk "BEGIN {exit !(${genome_ratio} < 0.75)}"; then 
        echo -n "WARN:genome ratio<0.75, " >> QC_EVAL
    fi

    # write pass if no fail
    if [ -f QC_EVAL ]; then
        # remove trailing comma
        sed -i '\$ s/, \$//; s/\\n//' QC_EVAL
    else
        echo "PASS" | tee QC_EVAL    
    fi

    qc_eval=\$(cat QC_EVAL)
    """
}


