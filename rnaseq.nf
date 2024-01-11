#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    Workflow parameters
========================================================================================
*/
params.input = "/vcu_gpfs2/home/mccbnfolab/balami/projects/wang_shawn/wang_shawn231101/raw_data/*_{R1,R2}.fastq.gz"
params.adaptor = "/vcu_gpfs2/home/morecockcm/bin/trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10"
params.genome = "/vcu_gpfs2/home/morecockcm/rnaseq_pipeline/genome_indexes/mouse/GRCm39_mm39/ref/"
params.gtf_file = "/vcu_gpfs2/home/balami/orfseq/reference/Streptococcus_sanguinis_SK36R_6908.current.gtf"
params.species = "MOUSE"
params.outdir = "ORFseq_results"

log.info """\
    R N A S E Q   P I P E L I N E
    ===================================
    reads                : ${params.input}
    output dir           : ${params.outdir}
    current adapter      : ${params.adaptor}
    annotation           : ${params.gtf_file}
    species              : ${params.species}
    """
    .stripIndent()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE PROCESSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


process TRIMMOMATIC {
    tag "Trimming unwanted bases"
    cpus 16
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2':
        'biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    input:
    path reads
    path adaptor

    output:
    path("*.paired.trim*.fastq.gz")   , emit: trimmed_reads
    path("*.unpaired.trim_*.fastq.gz"), optional:true, emit: unpaired_reads
    path("*.log")                     , emit: log
    path("*.summary")                 , emit: summary
    path "versions.yml"               , emit: versions

    script:
    def prefix = ${reads.baseName}

    """
    trimmomatic $trimmed -threads $task.cpus -trimlog ${prefix}.log \\
        -summary ${prefix}.summary \\
        ILLUMINACLIP:${adaptor}
        $reads \\
        $output \\
        $qual_trim \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version)
    END_VERSIONS
    """
}


process SUBREAD_FEATURECOUNTS {
    cpus 16

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0' :
        'biocontainers/subread:2.0.1--hed695b0_0' }"

    input:
    path bams
    path annotation

    output:
    path "*featureCounts.txt"         , emit: counts
    path "*featureCounts.txt.summary" , emit: summary
    path "versions.yml"               , emit: versions

    script:
    def prefix = ${bam.baseName}

    """
    featureCounts -T $task.cpus -F GTF -t exon -g gene_id -p --countReadPairs \\
        -a $annotation -o ${prefix}.featureCounts.txt \\
        ${bams.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g")
    END_VERSIONS
    """
}


process QUALIMAP_RNASEQ {
    tag "Qualimap BamQC"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qualimap:2.2.2d--1' :
        'biocontainers/qualimap:2.2.2d--1' }"

    input:
    path bam
    val species

    output:
    path "${prefix}"                  , emit: results
    path  "versions.yml"              , emit: versions

    script:
    def prefix = ${bam.baseName}

    """
    qualimap \\
        --java-mem-size=2G \\
        bamqc \\
        -bam $bam \\
        -gd $species \\
        -nt 12 -c \\
        -outfile ${prefix}_bamqc.pdf -outformat PDF:HTML \\
        -outdir $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
    END_VERSIONS
    """
}


process STAR_ALIGN {
    tag "Aligning using Star"
    cpus 16

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' :
        'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"

    input:
    path reads
    path index

    output:
    path '*sortedByCoord.out.bam'   , emit: bam_sorted
    path '*Log.final.out'           , emit: log_final
    path '*Log.out'                 , emit: log_out
    path '*Log.progress.out'        , emit: log_progress
    path  "versions.yml"            , emit: versions

    path '*d.out.bam'               , optional:true, emit: bam
    path '*fastq.gz'                , optional:true, emit: fastq
    path '*.out.sam'                , optional:true, emit: sam

    script:
    def prefix = ${reads.baseName}
    def out_sam_type = --outSAMtype BAM SortedByCoordinate
    
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn ${reads} \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        --readFilesCommand zcat \\
        $out_sam_type

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}Xd.out.bam
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.sortedByCoord.out.bam
    touch ${prefix}.toTranscriptome.out.bam
    touch ${prefix}.Aligned.unsort.out.bam
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.tab
    touch ${prefix}.SJ.out.tab
    touch ${prefix}.ReadsPerGene.out.tab
    touch ${prefix}.Chimeric.out.junction
    touch ${prefix}.out.sam
    touch ${prefix}.Signal.UniqueMultiple.str1.out.wig
    touch ${prefix}.Signal.UniqueMultiple.str1.out.bg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
}


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
    Channel
        .fromPath(params.genome, checkIfExists: true)
        .set {index_ch}
    TRIMMOMATIC(ch_fastq ,params.adaptor)
    trim_ch = TRIMMOMATIC.out.trimmed_reads
    STAR_ALIGN(trim_ch, index_ch)
    align_ch = STAR_ALIGN.out.bam_sorted
    featurecounts_ch = SUBREAD_FEATURECOUNTS(align_ch, params.gtf_file)
    qualimap_ch = QUALIMAP_RNASEQ(align_ch, params.species)

    // pre processing QC
    preTrimFastqc_ch = FASTQC(ch_fastq)
    MULTIQC((preTrimFastqc_ch).collect())

    // post trim QC
    postTrimFastqc_ch = FASTQC(trim_ch)
    MULTIQC((postTrimFastqc_ch).collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}