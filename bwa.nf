nextflow.enable.dsl=2

samples_ch = Channel.fromPath(params.samples, checkIfExists: true)
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
        bwa mem -t $task.cpus \
        ${ref_gen}/${params.ref_gen_prefix} \
        $fastq_1 $fastq_2 \
        2> ${sample_id}.err | \
        samtools view -@ $task.cpus -bT \
        ${ref_gen}/${params.ref_gen_prefix}.fa \
        > $bam
        """
}

process sort {
    publishDir "$params.outdir", pattern: "$bam_sorted", mode: "copy"

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path(bam_sorted)

    script:
        bam_sorted = "${sample_id}.sorted.bam"

        """
        samtools sort -@ $task.cpus -O bam -o $bam_sorted -T temp $bam 
        """
}

process bam_index {
    publishDir "$params.outdir", pattern: "$idx", mode: "copy"

    input:
        tuple val(sample_id), path(bam_sorted)

    output:
        tuple val(sample_id), path(bam_sorted), path(idx)

    script:
        idx = "${bam_sorted}.bai"

        """
        samtools index -@ $task.cpus $bam_sorted
        """
}

process flagstat {
    publishDir "$params.outdir", mode: "copy"

    input:
        tuple val(sample_id), path(sorted)
    
    output:
        path stat

    script:
        stat = "${sample_id}.flagstat.tsv"

        // There is also a more comprehensive samtools 'stats' command.
        """
        samtools flagstats -@ $task.cpus -O tsv $sorted > $stat
        """
}

workflow {
    samples_ch.combine(ref_gen_ch) | \
    bwa_mem | \
    sort | \
    bam_index | \
    map( it -> it[0..1] ) | \
    flagstat
}