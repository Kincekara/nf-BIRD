process PLASMIDFINDER {
    label 'process_low'
    container 'kincekara/plasmidfinder:2.1.6-db_2024-11-14'

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
    plasmidfinder.py \\
        -i ${assembly} \\
        -p /database/ \\
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
        plasmidfinder: \$(cat /VERSION)
        plasmidfinder_db: \$(cat /DB_DATE)
    END_VERSIONS
    """
}
