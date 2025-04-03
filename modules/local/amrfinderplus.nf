process AMRFINDER {
    label 'process_medium'
    container 'staphb/ncbi-amrfinderplus:4.0.15-2024-12-18.1'

    input:
    tuple val(samplename), path(scaffolds), path(prodigal_faa), path(prodigal_gff), val(bracken_organism), val(mash_organism)

    output:
    path "*.tsv"
    tuple val(samplename), env('amr_genes'), env('stress_genes'), env('virulence_genes'), env('amr_subclass'), emit: amr
    tuple val(samplename), path("*_amrfinder_all.tsv"), emit: report
    path "versions.yml", emit: versions

    script:
    def assembly = scaffolds
    """
    # select mash organism if avalible
    if [[ "${mash_organism}" != "null" ]]; then
        organism="${mash_organism}"
    else
        organism="${bracken_organism}"
    fi
    echo "organism is set to: \$organism"

    ## curated organisms ##
    # A. baumannii-calcoaceticus species complex
    declare -a abcc=(
        "Acinetobacter baumannii"
        "Acinetobacter calcoaceticus"
        "Acinetobacter lactucae"
        "Acinetobacter nosocomialis"
        "Acinetobacter pittii"
        "Acinetobacter seifertii"
    )
    # Burkholderia cepacia species complex
    declare -a bcc=(
        "Burkholderia aenigmatica"
        "Burkholderia ambifaria"   
        "Burkholderia anthina"   
        "Burkholderia arboris"   
        "Burkholderia catarinensis"   
        "Burkholderia cenocepacia"   
        "Burkholderia cepacia" 
        "Burkholderia cf. cepacia"  
        "Burkholderia contaminans"   
        "Burkholderia diffusa"   
        "Burkholderia dolosa"   
        "Burkholderia lata"   
        "Burkholderia latens"
        "Burkholderia metallica"  
        "Burkholderia multivorans"   
        "Burkholderia orbicola"   
        "Burkholderia paludis"   
        "Burkholderia pseudomultivorans"   
        "Burkholderia puraquae"   
        "Burkholderia pyrrocinia"   
        "Burkholderia semiarida"   
        "Burkholderia seminalis"   
        "Burkholderia sola"   
        "Burkholderia stabilis"   
        "Burkholderia stagnalis"   
        "Burkholderia territorii"   
        "Burkholderia ubonensis"   
        "Burkholderia vietnamiensis" 
    )
    # Burkholderia pseudomallei species complex
    declare -a bpc=(
        "Burkholderia humptydooensis"   
        "Burkholderia mallei"   
        "Burkholderia mayonis"   
        "Burkholderia oklahomensis"   
        "Burkholderia pseudomallei"   
        "Burkholderia savannae"   
        "Burkholderia singularis"   
        "Burkholderia thailandensis"   
    )
    # other species
    declare -a taxa=(   
        "Citrobacter freundii"
        "Clostridioides difficile"
        "Enterobacter asburiae"
        "Enterobacter cloacae"
        "Enterococcus faecalis"
        "Haemophilus influenzae"    
        "Klebsiella oxytoca"
        "Neisseria meningitidis"
        "Neisseria gonorrhoeae"
        "Pseudomonas aeruginosa" 
        "Serratia marcescens"  
        "Staphylococcus aureus"
        "Staphylococcus pseudintermedius"
        "Streptococcus agalactiae"
        "Streptococcus pyogenes"
        "Vibrio cholerae"
        "Vibrio parahaemolyticus"
        "Vibrio vulnificus"
    )

    # check organism in curated organism list
    genus=\$(echo \$organism | cut -d " " -f1)
    taxon=\$(echo \$organism | cut -d " " -f1,2)

    amrfinder_organism=""

    if [[ "\$genus" == "Acinetobacter" ]]; then
        for i in "\${abcc[@]}"; do
            if [[ "\$taxon" == "\$i" ]]; then
                amrfinder_organism="Acinetobacter_baumannii"
                break
            fi
        done
    elif [[ "\$genus" == "Burkholderia" ]]; then
        for i in "\${bcc[@]}"; do
            if [[ "\$taxon" == "\$i" ]]; then
                amrfinder_organism="Burkholderia_cepacia"
                break
            fi
        done
        for i in "\${bpc[@]}"; do
            if [[ "\$taxon" == "\$i" ]]; then
                amrfinder_organism="Burkholderia_pseudomallei"
                break
            fi
        done
    elif [[ "\$genus" == "Shigella" ]] || [[ "\$genus" == "Escherichia" ]]; then
        amrfinder_organism="Escherichia"
    elif [[ "\$genus" == "Salmonella" ]]; then
        amrfinder_organism="Salmonella"
    elif [[ "\$taxon" == "Campylobacter coli" ]] || [[ "\$taxon" == "Campylobacter jejuni" ]]; then
        amrfinder_organism="Campylobacter"
    elif [[ "\$taxon" == "Enterococcus faecium" ]] || [[ "\$taxon" == "Enterococcus hirae" ]]; then
        amrfinder_organism="Enterococcus_faecium"
    elif [[ "\$taxon" == "Klebsiella pneumoniae" ]] || [[ "\$taxon" == "Klebsiella aerogenes" ]]; then
        amrfinder_organism="Klebsiella_pneumoniae"
    elif [[ "\$taxon" == "Streptococcus pneumoniae" ]] || [[ "\$taxon" == "Streptococcus mitis" ]]; then
        amrfinder_organism="Streptococcus_pneumoniae"
    else    
        for i in "\${taxa[@]}"; do
            if [[ "\$taxon" == "\$i" ]]; then
                amrfinder_organism=\${taxon// /_}
                break
            fi
        done
    fi

    # checking bash variable
    echo "amrfinder_organism is set to: \$amrfinder_organism"
    
    # protein + nucleotide (activate HMM)
    if [[ -f "${prodigal_faa}" ]] && [[ -f "${prodigal_gff}" ]]; then
        # protein + nucleotide & use --organism 
        if [[ -n \$amrfinder_organism ]]; then
            amrfinder \\
                --plus \\
                --organism "\$amrfinder_organism" \\
                --name ${samplename} \\
                --nucleotide ${assembly} \\
                --protein ${prodigal_faa} \\
                --gff ${prodigal_gff} \\
                --annotation_format prodigal \\
                -o "${samplename}_amrfinder_all.tsv" \\
                --threads ${task.cpus} 2>&1 | tee amrfinder.STDOUT-and-STDERR.log
        # protein + nucleotide & no organism
        else        
            amrfinder \\
                --plus \\
                --name ${samplename} \\
                --nucleotide ${assembly} \\
                --protein ${prodigal_faa} \\
                --gff ${prodigal_gff} \\
                --annotation_format prodigal \\
                -o "${samplename}_amrfinder_all.tsv" \\
                --threads ${task.cpus} 2>&1 | tee amrfinder.STDOUT-and-STDERR.log
        fi
    # nucleotide only
    else
        # nucletode only & use --organism 
        if [[ -n \$amrfinder_organism ]] ; then
            amrfinder \\
            --plus \\
            --organism "\$amrfinder_organism" \\
            --name ${samplename} \\
            --nucleotide ${assembly} \\
            -o "${samplename}_amrfinder_all.tsv" \\
            --threads ${task.cpus} 2>&1 | tee amrfinder.STDOUT-and-STDERR.log
        # nucletode only & no organism 
        else 
            amrfinder \\
            --plus \\
            --name ${samplename} \\
            --nucleotide ${assembly} \\
            -o "${samplename}_amrfinder_all.tsv" \\
            --threads ${task.cpus} 2>&1 | tee amrfinder.STDOUT-and-STDERR.log
        fi
    fi

    # Element Type possibilities: AMR, STRESS, and VIRULENCE 
    # create headers for 3 output files; tee to 3 files and redirect STDOUT to dev null so it doesn't print to log file
    head -n 1 ${samplename}_amrfinder_all.tsv | tee ${samplename}_amrfinder_stress.tsv ${samplename}_amrfinder_virulence.tsv ${samplename}_amrfinder_amr.tsv >/dev/null
    # looks for all rows with STRESS, AMR, or VIRULENCE and append to TSVs
    # || true is so that the final grep exits with code 0, preventing failures
    grep 'STRESS' ${samplename}_amrfinder_all.tsv >> ${samplename}_amrfinder_stress.tsv || true
    grep 'VIRULENCE' ${samplename}_amrfinder_all.tsv >> ${samplename}_amrfinder_virulence.tsv || true
    grep 'AMR' ${samplename}_amrfinder_all.tsv >> ${samplename}_amrfinder_amr.tsv || true

    # create string outputs for all genes identified in AMR, STRESS, VIRULENCE
    amr_genes=\$(awk -F '\\t' '{ print \$7 }' ${samplename}_amrfinder_amr.tsv | tail -n+2 | tr '\\n' ', ' | sed 's/.\$//')
    stress_genes=\$(awk -F '\\t' '{ print \$7 }' ${samplename}_amrfinder_stress.tsv | tail -n+2 | tr '\\n' ', ' | sed 's/.\$//')
    virulence_genes=\$(awk -F '\\t' '{ print \$7 }' ${samplename}_amrfinder_virulence.tsv | tail -n+2 | tr '\\n' ', ' | sed 's/.\$//')
    amr_subclass=\$(awk -F '\\t' '{ print \$13 }' ${samplename}_amrfinder_amr.tsv | tail -n+2 | tr '\\n' ', ' | sed 's/.\$//')

    # if variable for list of genes is EMPTY, write string saying it is empty to float to Terra table
    if [[ -z "\$amr_genes" ]]; then
        amr_genes="No AMR genes detected by NCBI-AMRFinderPlus"
    fi 
    if [[ -z "\$stress_genes" ]]; then
        stress_genes="No STRESS genes detected by NCBI-AMRFinderPlus"
    fi 
    if [[ -z "\$virulence_genes" ]]; then
        virulence_genes="No VIRULENCE genes detected by NCBI-AMRFinderPlus"
    fi
    if [[ -z "\$amr_subclass" ]]; then
        amr_subclass="No AMR detected by NCBI-AMRFinderPlus"
    fi

    # version control
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus_db: \$(amrfinder --database_version | grep "Database version" | awk -F': ' '{print \$2}')
    END_VERSIONS
    """
}