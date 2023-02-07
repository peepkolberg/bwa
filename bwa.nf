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

process bam_index {
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

    input:
        tuple val(sample_id), path(sorted), path(idx)
    
    output:
        tuple path(sorted), path(idx), path(stat)

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