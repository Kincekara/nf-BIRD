process PROFILE {
    container 'staphb/bracken:2.9'

    input:
    tuple val(samplename), path(reads)
    path kraken2_db
    
    output:
    path "*.txt"
    tuple val(samplename), env('GENUS'), emit: genus
    tuple val(samplename), env('TAXON'), emit: taxon
    tuple val(samplename), env('TAXON'), env('RATIO'), emit: stats
    tuple val(samplename), path("*.bracken.filtered.txt"), path("*.taxid.txt"), emit: report
    path "versions.yml", emit: versions

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    def args = task.ext.args
    def args2 = task.ext.args2

    """
    # Decompress the Kraken2 database
    mkdir db
    tar -I pigz -C ./db/ -xvf ${kraken2_db} 

    # Run Kraken2
    kraken2 \\
    --db ./db/ \\
    --threads ${task.cpus} \\
    --report ${samplename}.kraken.report.txt \\
    --gzip-compressed \\
    --paired \\
    ${args} \\
    --report-minimizer-data \\
    ${read1} ${read2}
    
    # Run bracken
    bracken \\
    -d ./db/ \\
    -i ${samplename}.kraken.report.txt \\
    -o ${samplename}.bracken.txt \\
    -l S \\
    ${args2}

    # filter report
    awk 'NR==1; NR>1 {if (\$NF >= 0.01){print}}' ${samplename}.bracken.txt > ${samplename}.bracken.filtered.txt
    TAXON=\$(awk '{print \$NF,\$0}' ${samplename}.bracken.txt | sort -nr | cut -f2- -d' ' | awk -F "\\t" 'NR==1 {print \$1}')
    RATIO=\$(awk '{print \$NF,\$0}' ${samplename}.bracken.txt | sort -nr | cut -f2- -d' ' | awk -F "\\t" 'NR==1 {printf "%.2f\\n", \$NF*100}')
    awk '{print \$NF,\$0}' ${samplename}.bracken.txt | sort -nr | cut -f2- -d' ' | awk -F "\\t" 'NR==1 {print \$2}' > ${samplename}.taxid.txt 
    GENUS=\$(awk '{print \$NF,\$0}' ${samplename}.bracken.txt | sort -nr | cut -f2- -d' ' | awk 'NR==1 {print \$1}')

    # version control
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        kraken2_db: \$(basename -s .tar.gz ${kraken2_db} | cut -d "_" -f2,3,4)
        bracken: \$(echo \$(bracken -v) | cut -f2 -d'v')
    END_VERSIONS
    """

}

