process MLST {
    container 'staphb/mlst:2.25.0'

    input:
    tuple val(samplename), path(scaffolds)

    output:
    path "*_mlst.tsv"
    tuple val(samplename), path("*_mlst.tsv"), emit: report
    tuple val(samplename), env('predicted_mlst'), env('pubmlst_scheme'), emit: mlst
    path "versions.yml", emit: versions

    script:
    def assembly = scaffolds
    """  
    # create output header
    echo -e "Filename\\tPubMLST_Scheme_name\\tSequence_Type_(ST)\\tAllele_IDs" > ${samplename}_ts_mlst.tsv
    
    mlst \\
        --threads ${task.cpus} \\
        --nopath \\
        ${assembly} \\
        >> ${samplename}_ts_mlst.tsv
    
    # parse ts mlst tsv for relevant outputs
    if [ "\$(tail -n +2 ${samplename}_ts_mlst.tsv | wc -l)" -eq 0 ]; then
        predicted_mlst="No ST predicted"
        pubmlst_scheme="NA"
    else
        pubmlst_scheme="\$(cut -f2 ${samplename}_ts_mlst.tsv | tail -n 1)"
        predicted_mlst="ST\$(cut -f3 ${samplename}_ts_mlst.tsv | tail -n 1)"
        if [ "\$pubmlst_scheme" = "-" ]; then
            predicted_mlst="No ST predicted"
            pubmlst_scheme="NA"
        else
            if [ "\$predicted_mlst" = "ST-" ]; then
            predicted_mlst="No ST predicted"
            fi
        fi  
    fi
    
    # version control
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
    END_VERSIONS
    """
}