process PREDICT_TAXON {
    container 'staphb/mash:2.3-CBIRDv2'

    input:
    tuple val(samplename), path(scaffolds), val(genus)

    output:
    path "*.mash.sorted.tsv", emit: screen
    tuple val(samplename), path("*.top_taxon.tsv"), emit: report
    tuple val(samplename), env('TAXON'), emit: taxon
    tuple val(samplename), env('TAXON'), env('RATIO'), emit:stats
    path "versions.yml", emit: versions

    script:
    def assembly = scaffolds
    """
    # screen assembly
    mash screen -p ${task.cpus} /db/cbird-v2.0-lite.msh ${assembly} > ${samplename}.mash.tsv

    # parse results
    sort -gr ${samplename}.mash.tsv > ${samplename}.mash.sorted.tsv
    top=\$(awk -F "\\t" 'NR==1 {print \$6}' ${samplename}.mash.sorted.tsv)
    if [[ "\$top" =~ ^\\[[0-9]+[[:space:]]seqs\\]+ ]]; then
        candidate=\$(echo "\$top" | cut -d ' ' -f 4,5,6,7)
    else
        candidate=\$(echo "\$top" | cut -d ' ' -f 2,3,4,5)
    fi
    # check subspecies
    if [[ \$(echo "\$candidate" | awk '{print \$3}') == "subsp." ]]; then
        TAXON="\$candidate"
    else
        TAXON=\$(echo "\$candidate" | cut -d ' ' -f 1,2)
    fi

    # find ratio
    RATIO=\$(awk -F "\\t" 'NR==1 {printf "%.2f\\n",\$1*100}' ${samplename}.mash.sorted.tsv)
    printf '%s\\t%s\\n' "\$TAXON" "\$RATIO" > ${samplename}.top_taxon.tsv

    # version control
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$( mash --version )
    END_VERSIONS
    """
}