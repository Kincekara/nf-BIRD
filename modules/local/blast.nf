process TBLASTN {
    label 'process_low'
    container 'staphb/blast:2.15.0'

    input:
    tuple val(samplename), path(scaffolds)
    path target_genes_fasta

    output:
    path "*.tsv"
    tuple val(samplename), path("*.tblastn.tsv"), emit: report
    tuple val(samplename), env('GENES'), emit: stats
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args
    def args2 = task.ext.args2
    def query = target_genes_fasta
    def subject = scaffolds
    """
    if [ ! -f "${query}" ]
    then 
        echo "No query file!"
    elif [ ! -f "${subject}" ] 
    then
        echo "No subject file!"
    else
        tblastn -query ${query} -subject ${subject} ${args} -out blast.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
        awk 'BEGIN{print "qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen"};\$3>=${args2} {print}' blast.txt > ${samplename}.tblastn.tsv
        hits=\$(awk '\$3>=${args2} {print \$1}' blast.txt | tr '\\n' ', ' | sed 's/.\$//')
    fi
    if [ -z "\$hits" ]
    then
        echo "No blast hit!" > GENES
    else
        echo "\$hits" > GENES
    fi

    # version control
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(tblastn -version 2>&1 | sed 's/^.*tblastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
