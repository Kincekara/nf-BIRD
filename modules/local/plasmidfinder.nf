process PLASMIDFINDER {
    label 'process_low'
    container 'staphb/plasmidfinder:3.0.3'

    input:
    tuple val(samplename), path(scaffolds)

    output:
    path "*.tsv"
    tuple val(samplename), env('PF'), emit: plasmids
    tuple val(samplename), path("*.plasmid.tsv"), emit: report
    path "versions.yml", emit: versions

    script:
    def assembly = scaffolds
    def min_coverage = 0.6
    def threshold = 0.9
    """   
    # Run plasmidfinder
    python -m plasmidfinder \\
        -i ${assembly} \\
        -l ${min_coverage} \\
        -t ${threshold} \\
        -x

    # parse outputs
    if [ ! -f results_tab.tsv ]; then
        PF="No plasmids detected in database"
    else
        PF="\$(tail -n +2 results_tab.tsv | uniq | cut -f 2 | sort | paste -s -d, - )"
        if [ "\$PF" == "" ]; then
            PF="No plasmids detected in database"
        fi  
    fi
    echo "\$PF" | tee PLASMIDS

    # rename results
    mv results_tab.tsv ${samplename}.plasmid.tsv
    
    # version control
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plasmidfinder: \$(python -m plasmidfinder -v)
        plasmidfinder_db: \$(</plasmidfinder_db/VERSION)
    END_VERSIONS
    """
}
