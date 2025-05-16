process SPADES {
    label 'process_high_memory'
    container 'staphb/spades:4.2.0'

    input:
    tuple val(samplename), path(reads), val(read_length)

    output:
    tuple val(samplename), file("*_contigs.fasta"), emit: contigs
    tuple val(samplename), file("*_scaffolds_trim.fasta"), emit: scaffolds
    path "versions.yml", emit: versions

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    def contig_threshold = read_length.toInteger() * 2    
    """
    # version control
    spades.py -v > VERSION 

    # assembly
    spades.py \\
        --careful \\
        --only-assembler \\
        --pe1-1 ${read1} \\
        --pe1-2 ${read2} \\
        -o out

    # get & rename output   
    mv out/contigs.fasta ${samplename}_contigs.fasta
    mv out/scaffolds.fasta ${samplename}_scaffolds.fasta

    # remove short scaffolds
    echo "Removing scaffolds shorter than ${contig_threshold}"
    trim_scaffolds.py ${samplename}_scaffolds.fasta ${samplename}_scaffolds_trim.fasta ${contig_threshold}

    # version control
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed -n 's/^.*SPAdes genome assembler v//p')
    END_VERSIONS
    """



}