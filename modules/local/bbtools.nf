process ASSEMBLY_PREP {
    label 'process_medium'
    container 'staphb/bbtools:39.49'

    input:
    tuple val(samplename), path(reads_trimmed), val(total_reads) // from ch_filtered_reads

    output:
    tuple val(samplename), file("*_{1,2}.clean.norm.fastq.gz"), emit: reads
    tuple val(samplename), env('PHIX_RATIO'), emit: phix
    path "*.phix.stats.txt"
    path "*.bbnorm.log", optional: true
    path "versions.yml", emit: versions
    

    script:
    def read1 = reads_trimmed[0]
    def read2 = reads_trimmed[1]
    def normalization = task.ext.args.normalization
    def norm_target = task.ext.args.norm_target
    def min_depth = task.ext.args.min_depth
    def read_threshold = task.ext.args.read_threshold
    def memory = "-Xmx${Math.round(Math.max(1, Math.floor(task.memory.toGiga() * 0.95)))}g"
    """
    # PhiX cleaning   
    bbduk.sh \\
        in1=${read1} \\
        in2=${read2} \\
        out1=${samplename}_1.clean.fastq.gz \\
        out2=${samplename}_2.clean.fastq.gz \\
        outm=${samplename}.matched_phix.fq \\
        ref=/bbmap/resources/phix174_ill.ref.fa.gz \\
        k=31 \\
        hdist=1 \\
        bgzip=f \\
        unbgzip=f \\
        stats=${samplename}.phix.stats.txt

    PHIX_RATIO=\$(grep Matched ${samplename}.phix.stats.txt | awk '{print \$3}')
    
    # normalization
    if ${normalization} && [ "${total_reads}" -gt "${read_threshold}" ]      
    then
        echo "normalizing reads..."
        bbnorm.sh \\
            ${memory} \\
            threads=${task.cpus} \\
            in=${samplename}_1.clean.fastq.gz \\
            in2=${samplename}_2.clean.fastq.gz \\
            out=${samplename}_1.clean.norm.fastq.gz \\
            out2=${samplename}_2.clean.norm.fastq.gz \\
            target=${norm_target} \\
            min=${min_depth} \\
            &> ${samplename}.bbnorm.log
    else
        echo "skipping normalization..."
        mv ${samplename}_1.clean.fastq.gz ${samplename}_1.clean.norm.fastq.gz
        mv ${samplename}_2.clean.fastq.gz ${samplename}_2.clean.norm.fastq.gz
    fi

    # version control
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbtools: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
    
    
}





