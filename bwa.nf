nextflow.enable.dsl=2

samples_ch = Channel.fromPath(params.samples)
    .ifEmpty { error "Cannot find samples file: ${params.samples}" }
    .splitCsv(header: true, sep: "\t", strip: true)
    .map{row -> [row.sample_id, row.fastq_1, row.fastq_2]}

ref_gen_ch = Channel.fromPath(params.ref_gen, type: "dir")

process bwa_mem {
    input:
        tuple val(sample_id), path(fastq_1), path(fastq_2), path(ref_gen)

    output:
        tuple val(sample_id), path(bam)

    script:
        bam = "${sample_id}.bam"

        """
        bwa mem -t 8 \
        ${ref_gen}/${params.ref_gen_prefix} \
        $fastq_1 $fastq_2 \
        2> ${sample_id}.err | \
        samtools view -@ 8 -bT \
        ${ref_gen}/${params.ref_gen_prefix}.fa \
        > $bam
        """
}

process sort {
    publishDir "${projectDir}/results", pattern: "$bam_sorted", mode: "copy"

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path(bam_sorted)

    script:
        bam_sorted = "${sample_id}.sorted.bam"

        """
        samtools sort -@ 4 -O bam -o $bam_sorted -T temp $bam 
        """
}

// Index gets published by move (not copy) but that's ok because the next process doesn't use index anyway.
process bam_index {
    publishDir "${projectDir}/results", pattern: "$idx", mode: "move"

    input:
        tuple val(sample_id), path(bam_sorted)

    output:
        tuple val(sample_id), path(bam_sorted), path(idx)

    script:
        idx = "${bam_sorted}.bai"

        """
        samtools index -@ 4 $bam_sorted
        """
}

process flagstat {
    publishDir "${projectDir}/results", mode: "move"

    // Even though this process doesn't use the index file we created in the previous process,
    // it must still be declared in the input so the dimensionality of this input matches the
    // output of the previous process.
    input:
        tuple val(sample_id), path(sorted), path(idx)
    
    output:
        path stat

    script:
        stat = "${sample_id}.flagstat.tsv"

        // There is also a more comprehensive samtools 'stats' command.
        """
        samtools flagstats -@ 4 -O tsv $sorted > $stat
        """
}

workflow {
    samples_ch.combine(ref_gen_ch) | bwa_mem | sort | bam_index | flagstat
}