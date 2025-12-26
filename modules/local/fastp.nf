process FASTP {
    label 'process_medium'
    container 'staphb/fastp:1.0.1'

    input:
    tuple val(samplename), path(reads)
    path adapters

    output:
    tuple val(samplename), file("*_R{1,2}_trim.fastq.gz"), env('TOTAL_READS'), emit: reads
    path "*_fastp.html"
    path "*_fastp.json", emit: json
    tuple val(samplename), path ("*_total_bases.txt"), path("*_fastp.html"), emit: report
    tuple val(samplename), env('READ_LENGTH'), emit: read_length
    tuple val(samplename), env('R1_READS'), env('R2_READS'), env('TOTAL_READS'), env('TOTAL_READS_TRIM'), env('R1_Q30_RAW'), env('R1_Q30_TRIM'), env('R2_Q30_RAW'), env('R2_Q30_TRIM'), emit: stats
    tuple val(samplename), env('TOTAL_READS_TRIM'), env('R1_Q30_TRIM'), env('R2_Q30_TRIM'), emit: stats2
    path "versions.yml", emit: versions

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    def adapter_arg = (adapters && adapters.size() > 0)? "--adapter_fasta ${adapters}" : "--detect_adapter_for_pe"
    def args = task.ext.args.join(" ")
    """
    fastp \\
        -i ${read1} \\
        -I ${read2} \\
        -o ${samplename}_R1_trim.fastq.gz \\
        -O ${samplename}_R2_trim.fastq.gz \\
        ${args} \\
        ${adapter_arg} \\
        --thread ${task.cpus} \\
        -h ${samplename}_fastp.html

    # parse output
    jq '.summary.before_filtering.total_bases' fastp.json > ${samplename}_total_bases.txt
    jq '.summary.after_filtering.total_bases' fastp.json >> ${samplename}_total_bases.txt  
    R1_READS=\$(jq '.read1_before_filtering.total_reads' fastp.json)
    R2_READS=\$(jq '.read2_before_filtering.total_reads' fastp.json)
    TOTAL_READS=\$(jq '.summary.before_filtering.total_reads' fastp.json)
    TOTAL_READS_TRIM=\$(jq '.summary.after_filtering.total_reads' fastp.json)

    r1_q30=\$(jq '.read1_before_filtering.q30_bases' fastp.json)
    r1_total=\$(jq '.read1_before_filtering.total_bases' fastp.json)
    r2_q30=\$(jq '.read2_before_filtering.q30_bases' fastp.json)
    r2_total=\$(jq '.read2_before_filtering.total_bases' fastp.json)
    r1_q30_trim=\$(jq '.read1_after_filtering.q30_bases' fastp.json)
    r1_total_trim=\$(jq '.read1_after_filtering.total_bases' fastp.json)
    r2_q30_trim=\$(jq '.read2_after_filtering.q30_bases' fastp.json)
    r2_total_trim=\$(jq '.read2_after_filtering.total_bases' fastp.json)

    READ_LENGTH=\$((\$(jq '.summary.before_filtering.read1_mean_length' fastp.json) + 1 ))

    mv fastp.json ${samplename}_fastp.json
        
    if [ "\$r1_total_trim" -gt 0 ]
    then
        R1_Q30_RAW=\$(echo "\$r1_q30 \$r1_total" | awk '{printf "%.2f", (\$1/\$2)*100 }')
        R1_Q30_TRIM=\$(echo "\$r1_q30_trim \$r1_total_trim" | awk '{printf "%.2f", (\$1/\$2)*100 }')
    else
        R1_Q30_RAW=0
        R1_Q30_TRIM=0
    fi

    if [ "\$r2_total_trim" -gt 0 ]
    then
        R2_Q30_RAW=\$(echo "\$r2_q30 \$r2_total" | awk '{printf "%.2f", (\$1/\$2)*100 }')
        R2_Q30_TRIM=\$(echo "\$r2_q30_trim \$r2_total_trim" | awk '{printf "%.2f", (\$1/\$2)*100 }')
    else
        R2_Q30_RAW=0
        R2_Q30_TRIM=0
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}

