process MLST {
    container 'staphb/mlst:2.32.2'

    input:
    tuple val(samplename), path(scaffolds)

    output:
    path "*.mlst.tsv"
    tuple val(samplename), path("${samplename}.mlst.tsv"), emit: report
    tuple val(samplename), path("${samplename}.legacy.mlst.tsv"), emit: legacy_report
    tuple val(samplename), env('predicted_mlst'), env('pubmlst_scheme'), emit: mlst
    path "versions.yml", emit: versions

    script:
    def assembly = scaffolds
    """  
    # run mlst
    mlst \\
        --full \\
        --threads ${task.cpus} \\
        --nopath \\
        ${assembly} \\
        --outfile ${samplename}.mlst.tsv
    
    # rerun for A.baumannii oxford scheme
    scheme=\$(awk -F'\\t' 'NR==2 {print \$2}' ${samplename}.mlst.tsv)
    if [ "\$scheme" == "abaumannii_2" ]; then
        mlst \\
            --full \\
            --threads ${task.cpus} \\
            --nopath \\
            --scheme "abaumannii" \\
            ${assembly} \\
            --outfile ${samplename}.mlst.oxford.tsv
        # combine results
        awk 'NR==2' ${samplename}.mlst.oxford.tsv >> ${samplename}.mlst.tsv
    fi

    # parse ts mlst tsv for relevant outputs
    predicted_mlst=\$(awk -F'\\t' 'NR==2 {if (\$3 == "-") print "No ST predicted"; else print "ST"\$3 }' ${samplename}.mlst.tsv)
    pubmlst_scheme=\$(awk -F'\\t' 'NR==2 {if (\$2 == "-") print "NA"; else print \$2 }' ${samplename}.mlst.tsv)

    # create legacy output file for report scripts.
    echo -e "Filename\\tPubMLST_Scheme_name\\tSequence_Type_(ST)\\tAllele_IDs" > ${samplename}.legacy.mlst.tsv
    awk -F'\\t' 'NR==2 {print \$1"\\t"\$2"\\t"\$3"\\t"\$6}' ${samplename}.mlst.tsv >> ${samplename}.legacy.mlst.tsv

    
    # version control
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
    END_VERSIONS
    """
}