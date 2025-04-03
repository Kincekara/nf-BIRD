process QUAST {
    label 'process_single'
    container 'staphb/quast:5.3.0-slim'

    input:
    tuple val(samplename), path(scaffolds)

    output:
    path "*.tsv"
    tuple val(samplename), env('GENOME_LENGTH'), path("${samplename}_report.tsv"), emit: report
    tuple val(samplename), env('GENOME_LENGTH'), env('NUM_CONTIGS'), env('N50'), env('L90'), env('GC'), emit: stats
    tuple val(samplename), env('NUM_CONTIGS'), emit: num_contigs
    path "versions.yml", emit: versions

    script:
    def assembly = scaffolds
    """
    quast.py ${assembly} -o .   
    
    # parse results
    GENOME_LENGTH=\$(awk -F'\\t' '/^Total length\\t/ {print \$2}' report.tsv)
    NUM_CONTIGS=\$(awk -F'\\t' '/^# contigs\\t/ {print \$2}' report.tsv)
    N50=\$(awk -F'\\t' '/^N50\\t/ {print \$2}' report.tsv)
    L90=\$(awk -F'\\t' '/^L90\\t/ {print \$2}' report.tsv)
    GC=\$(awk -F'\\t' '/^GC \\(%\\)\\t/ {print \$2}' report.tsv)

    mv report.tsv ${samplename}_report.tsv

    # version control
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//; s/,//')
    END_VERSIONS
    """
}
