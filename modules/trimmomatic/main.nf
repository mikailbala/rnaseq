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
    trimmomatic \\
        $trimmed \\
        -threads $task.cpus \\
        -trimlog ${prefix}.log \\
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