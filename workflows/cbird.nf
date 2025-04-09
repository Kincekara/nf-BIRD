/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTP                  } from '../modules/local/fastp'
include { ASSEMBLY_PREP          } from '../modules/local/bbtools'
include { PROFILE                } from '../modules/local/taxonomy'
include { SPADES                 } from '../modules/local/spades'
include { QUAST                  } from '../modules/local/quast'
include { MLST                   } from '../modules/local/mlst'
include { AMRFINDER              } from '../modules/local/amrfinderplus'
include { PLASMIDFINDER          } from '../modules/local/plasmidfinder'
include { CHECKM2                } from '../modules/local/checkm2'
include { PREDICT_TAXON          } from '../modules/local/mash'
include { TBLASTN                } from '../modules/local/blast'
include { GENERATE_REPORT        } from '../modules/local/report'
include { QC_CHECK               } from '../modules/local/qc'
include { TABULATE               } from '../subworkflows/local/tabulate'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CBIRD {
    take:
    ch_samplesheet

    main:
    ch_versions = Channel.empty()

    // QC check, adapter removal, quality filtering and trimming
    FASTP(ch_samplesheet, params.adapters)
    ch_versions = ch_versions.mix(FASTP.out.versions)

    FASTP.out.reads
        .branch { it ->
            keep: it[2].toInteger() > params.minimum_total_reads
            remove: it[2].toInteger() <= params.minimum_total_reads
        }
        .set { reads_filter }
    reads_filter.remove.view { "${it[0]} : NOT ENOUGH READS, NOT BE ANALYZED!" }

    // Phix cleaning && normalization (if necessary)
    ASSEMBLY_PREP(reads_filter.keep)
    ch_versions = ch_versions.mix(ASSEMBLY_PREP.out.versions)

    // Taxonomic profiling and abundance estimation of reads
    PROFILE(ASSEMBLY_PREP.out.reads, params.kraken2_db)
    ch_versions = ch_versions.mix(PROFILE.out.versions)

    // De novo assembly
    SPADES(ASSEMBLY_PREP.out.reads.join(FASTP.out.read_length))
    ch_versions = ch_versions.mix(SPADES.out.versions)

    // Target gene search
    if (params.target_genes_fasta) {
        TBLASTN(SPADES.out.scaffolds, params.target_genes_fasta)
        ch_versions = ch_versions.mix(TBLASTN.out.versions)
    }

    // Species prediction for target genera
    def target_genera = [
        "Acinetobacter",
        "Burkholderia",
        "Citrobacter",
        "Enterobacter",
        "Escherichia",
        "Klebsiella",
        "Kluyvera",
        "Metapseudomonas",
        "Morganella",
        "Neisseria",
        "Proteus",
        "Providencia",
        "Pseudomonas",
        "Raoultella",
        "Salmonella",
        "Serratia",
        "Streptococcus",
    ]

    PROFILE.out.genus
        .branch { it ->
            cont: it[1] in target_genera
            skip: !(it[1] in target_genera)
        }
        .set { ch_genusfilter }

    ch_skipped = ch_genusfilter.skip.map { sample -> tuple(sample[0], 'NA', 'NA') }

    // Bacterial identification
    PREDICT_TAXON(SPADES.out.scaffolds.join(ch_genusfilter.cont))
    ch_versions = ch_versions.mix(PREDICT_TAXON.out.versions)

    // Genome assembly evaluation
    QUAST(SPADES.out.scaffolds)
    ch_versions = ch_versions.mix(QUAST.out.versions)

    // Completeness and contamination
    CHECKM2(SPADES.out.scaffolds, params.checkm2_db)
    ch_versions = ch_versions.mix(CHECKM2.out.versions)

    // MLST typing
    MLST(SPADES.out.scaffolds)
    ch_versions = ch_versions.mix(MLST.out.versions)

    // AMR gene detection
    AMRFINDER(
        SPADES.out.scaffolds
        .join(CHECKM2.out.prodigal)
        .join(PROFILE.out.taxon)
        .join(PREDICT_TAXON.out.taxon, remainder: true)
    )
    ch_versions = ch_versions.mix(AMRFINDER.out.versions)

    // Plasmid detection
    PLASMIDFINDER(SPADES.out.scaffolds)
    ch_versions = ch_versions.mix(PLASMIDFINDER.out.versions)

    // Individual summary report generation
    GENERATE_REPORT(
        FASTP.out.report
        .join(PROFILE.out.report)
        .join(MLST.out.report)
        .join(AMRFINDER.out.report)
        .join(PLASMIDFINDER.out.report)
        .join(ASSEMBLY_PREP.out.phix)
        .join(QUAST.out.report)
        .join(CHECKM2.out.report)
        .join(PREDICT_TAXON.out.report, remainder: true)
        .join(params.target_genes_fasta ? TBLASTN.out.report : Channel.empty(), remainder: true)
    )

    // Check QC metrics, mark passed and failed samples
    QC_CHECK(
        GENERATE_REPORT.out.stats
        .join(CHECKM2.out.stats)
        .join(QUAST.out.num_contigs)
        .join(FASTP.out.stats2)
    )

    // Aggragated summary report generation
    TABULATE(
        FASTP.out.stats,
        ASSEMBLY_PREP.out.phix,
        PROFILE.out.stats,
        QUAST.out.stats,
        CHECKM2.out.stats,
        PREDICT_TAXON.out.stats.mix(ch_skipped),
        MLST.out.mlst,
        AMRFINDER.out.amr,
        PLASMIDFINDER.out.plasmids,
        params.target_genes_fasta ? TBLASTN.out.stats : Channel.empty(),
        GENERATE_REPORT.out.stats,
        QC_CHECK.out.stats,
    )

    // Version
    softwareVersionsToYAML(ch_versions).collectFile(
        storeDir: "${params.outdir}/pipeline_info",
        name: "nfBIRD_v${workflow.manifest.version}" + '_versions.yml',
        sort: true,
        newLine: true,
    )

    emit:
    versions = ch_versions
}
