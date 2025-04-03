process CHECKM2 {
    label 'process_high'
    container 'staphb/checkm2:1.0.2'

    input:
    tuple val(samplename), path(scaffolds)
    path checkm2_db

    output:
    path "*.tsv"
    tuple val(samplename), path("*.report.tsv"), emit: report
    tuple val(samplename), path("*.faa"), path("*.gff"), emit: prodigal, optional: true
    tuple val(samplename), env('COMPLETENESS'), env('CONTAMINATION'), emit: stats
    path "versions.yml", emit: versions

    script:
    def assembly = scaffolds
    """
    # prep inputs
    tar -C . -xvf ${checkm2_db}     
    mkdir bins
    cp ${assembly} ./bins/
    
    # run chekm2    
    export TMPDIR=/tmp # -> https://github.com/broadinstitute/cromwell/issues/3647    
    
    checkm2 predict \\
        --threads ${task.cpus} \\
        --tmpdir /tmp \\
        -x fasta \\
        --input ./bins \\
        --output-directory ./out \\
        --database_path ./CheckM2_database/uniref100.KO.1.dmnd
    
    # parse results
    cp ./out/quality_report.tsv ./${samplename}.checkm2.report.tsv
    COMPLETENESS=\$(awk -F '\\t' 'NR==2 { print \$2 }' out/quality_report.tsv)
    CONTAMINATION=\$(awk -F '\\t' 'NR==2 { print \$3 }' out/quality_report.tsv)
        
    # (bonus) run prodigal to get protein data
    prodigal -m -i ${assembly} -f gff -o ${samplename}.prodigal.gff -a ${samplename}.prodigal.faa

    # version control
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """
}