process PLASMIDFINDER {
    label 'process_low'
    container 'staphb/plasmidfinder:3.0.1'

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
        -j data.json 

    # Create legacy tsv output
    json_to_tsv.py
    PF=\$(<PLASMIDS)

    # rename results
    mv results_tab.tsv ${samplename}.plasmid.tsv
    
    # version control
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plasmidfinder: \$(python -m plasmidfinder -v)
        plasmidfinder_db: \$(</database/VERSION.txt)
    END_VERSIONS
    """
}
