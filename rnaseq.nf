#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    Workflow parameters
========================================================================================
*/
params.input = "/vcu_gpfs2/home/balami/orfseq/data/*_{1,2}.fastq.gz"
params.adaptor = "/vcu_gpfs2/home/morecockcm/bin/trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10"
params.genome = "/vcu_gpfs2/home/morecockcm/rnaseq_pipeline/genome_indexes/mouse/GRCm39_mm39/ref/"
params.gtf_file = "/vcu_gpfs2/home/balami/orfseq/reference/Streptococcus_sanguinis_SK36R_6908.current.gtf"
params.species = "MOUSE"
params.outdir = "ORFseq_results"

log.info """\
    R N A S E Q   P I P E L I N E
    ===================================
    reads                : ${params.reads}
    output dir           : ${params.outdir}
    current adapter      : ${params.adaptor}
    annotation           : ${params.gtf_file}
    """
    .stripIndent()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { TRIMMOMATIC                 } from './modules/trimmomatic'
include { STAR_ALIGN                  } from './modules/star'
include { SUBREAD_FEATURECOUNTS       } from './modules/subread/featurecounts'
inlcude { QUALIMAP_RNASEQ             } from './modules/qualimap'


process FASTQC {
    tag "FASTQC on samples"
    cpus 12
    publishDir "$params.outdir", mode:'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    path reads

    output:
    path "FastQC"

    script:
    """
    mkdir FastQC
    fastqc --threads $task.cpus -o FastQC -q ${reads}
    """
}

process MULTIQC {
    tag "Producing a MultiQC Report"
    publishDir "$params.outdir", mode:'copy'
    cpus 12

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.17--pyhdfd78af_0' :
        'biocontainers/multiqc:1.17--pyhdfd78af_0' }"

    input:
    path '*.fastq'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {
    //
    // Create input channel from input file provided through params.input
    //
    Channel
        .fromFilePairs(params.input, checkIfExists: true)
        .set { ch_fastq }
    
    index_ch = fromPath(param.genome)
    trim_ch = TRIMMOMATIC(ch_fastq ,params.adaptor)
    align_ch = STAR_ALIGN(trim_ch.out[0], index_ch)
    featurecounts_ch = SUBREAD_FEATURECOUNTS(align_ch.out[0], parmas.gtf_file)
    qualimap_ch = QUALIMAP_RNASEQ(align_ch.out[0], params.species)

    // pre processing QC
    preTrimFastqc_ch = FASTQC(ch_fastq)
    MULTIQC((preTrimFastqc_ch).collect())

    // [post trim QC
    postTrimFastqc_ch = FASTQC(trim_ch)
    MULTIQC((postTrimFastqc_ch).collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}