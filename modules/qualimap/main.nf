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
    prefix   = ${bam.baseName}

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