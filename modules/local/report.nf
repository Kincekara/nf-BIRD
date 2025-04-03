process GENERATE_REPORT {
    label 'process_single'
    container 'staphb/cbird-util:2.0'

    input: 
    tuple val(samplename), 
        path (total_bases), 
        path(fastp_report), 
        path(taxon_report), 
        path(taxid), 
        path(mlst_report), 
        path(amr_report), 
        path(plasmid_report), 
        val(phix_ratio), 
        val(genome_length), 
        path(quast_report), 
        path(checkm2_report), 
        val(mash), 
        val(blast)


    output:
    path "*.html"
    path "*.docx", optional: true
    tuple val(samplename), env('genome_ratio'), env('coverage'), env('coverage_trim'), emit: stats

    script:
    def labid = samplename
    def args = task.ext.args.join(" ")
    def footer_note = task.ext.args2
    def analysis_date = task.ext.args3.date
    def version = task.ext.args3.version
    def mash_result = mash? "${mash}" : ""
    def blast_result = blast? "${blast}" : ""
    """
    # check mash results  
    if [ -f "${mash_result}" ]
    then
        # catch taxon & find genome size
        taxon=\$(awk -F '\\t' '{print \$1}' ${mash_result})
        percent=\$(awk -F '\\t' '{print \$2}' ${mash_result})
        datasets summary genome taxon "\$taxon" --reference > gs.json

        #  create plain report
        if [ -n "${labid}" ]
        then
            plain_report.py \\
            -d "${analysis_date}" \\
            -i "${labid}" \\
            -o "\$taxon" \\
            -p "\$percent" \\
            -a ${amr_report} \\
            ${args}
        else
            echo "No labid is provided. Skipping plain report generation."
        fi

        # create summary report with mash
        if [ -f "${blast_result}" ]
        # with blast
        then
            html_report.py \\
            -s ${samplename} \\
            -t ${taxon_report} \\
            -st ${mlst_report} \\
            -a ${amr_report} \\
            -p ${plasmid_report} \\
            -c "${version}" \\
            -m ${mash_result} \\
            -b ${blast_result} \\
            -f "${footer_note}"
        # mash only
        else
            html_report.py \\
            -s ${samplename} \\
            -t ${taxon_report} \\
            -st ${mlst_report} \\
            -a ${amr_report} \\
            -p ${plasmid_report} \\
            -c "${version}" \\
            -m ${mash_result} \\
            -f "${footer_note}"
        fi    
    else
    # find genome size with bracken taxon id
    taxid=\$(cat ${taxid})
    datasets summary genome taxon "\$taxid" --reference > gs.json
    
        # create summary report w/o mash      
        if [ -f "${blast_result}" ]
        # blast only
        then
            html_report.py \\
            -s ${samplename} \\
            -t ${taxon_report} \\
            -st ${mlst_report} \\
            -a ${amr_report} \\
            -p ${plasmid_report} \\
            -c "${version}" \\
            -b ${blast_result} \\
            -f "${footer_note}"
        # no mash or blast
        else
            html_report.py \\
            -s ${samplename} \\
            -t ${taxon_report} \\
            -st ${mlst_report} \\
            -a ${amr_report} \\
            -p ${plasmid_report} \\
            -c "${version}" \\
            -f "${footer_note}"
        fi
    fi

    # alternative source for expected genome size
    jq -r '.reports[0].assembly_stats.total_sequence_length' gs.json > alt_gs.txt
    
    # calculate esimated coverage & genome ratio
    taxid=\$(cat ${taxid})
    est_coverage.py \\
    ${total_bases} \\
    "\$taxid" \\
    alt_gs.txt \\
    "${genome_length}"    

    # create QC summary
    qc_report.py \\
    ${samplename} \\
    ${fastp_report} \\
    ${taxon_report} \\
    ${quast_report} \\
    ${checkm2_report} \\
    "${version}" \\
    "${phix_ratio}" \\
    "COVERAGE" \\
    "GENOME_RATIO"

    # Collect stats
    genome_ratio=\$(cat GENOME_RATIO)
    coverage=\$(cat COVERAGE)
    coverage_trim=\$(cat COVERAGE_TRIM)

    """
}
