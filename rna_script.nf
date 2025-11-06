#!/usr/bin/env nextflow

params.input_dir   = "${params.input_dir}"
params.output_dir  = "${params.output_dir}"
params.genomeIndex = "${params.genomeIndex}"
params.annotation  = "${params.annotation}"
params.adapters    = "${params.adapters}"
params.threads     = "${params.threads}"

process trimming {
    input:
        path reads
    output:
        path "${reads.simpleName}_trimmed.fastq.gz"
    script:
        """
        OUTPUT="${reads.simpleName}_trimmed.fastq.gz"
        trimmomatic SE -phred33 ${reads} $OUTPUT \
            ILLUMINACLIP:${params.adapters}:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """
}

process fastqc_trimmed {
    input:
        path trimmed_reads
    output:
        path "${trimmed_reads.simpleName}_fastqc.zip"
        path "${trimmed_reads.simpleName}_fastqc.html"
    script:
        """
        fastqc ${trimmed_reads} --outdir .
        """
}

process alignment {
    input:
        path trimmed_reads
    output:
        path "${trimmed_reads.simpleName}_Aligned.sortedByCoord.out.bam"
    script:
        """
        STAR \
            --genomeDir ${params.genomeIndex} \
            --readFilesIn ${trimmed_reads} \
            --runThreadN ${params.threads} \
            --outFileNamePrefix ${trimmed_reads.simpleName}_ \
            --outSAMtype BAM SortedByCoordinate
        """
}

process quantification {
    input:
        path bam_files
    output:
        path "gene_counts.txt"
    script:
        """
        featureCounts -a ${params.annotation} -o gene_counts.txt ${bam_files} -p -t exon -g gene_id -T ${params.threads}
        """
}

process multiqc {
    input:
        path all_results
    output:
        path "multiqc_report.html"
    script:
        """
        multiqc . -f -o .
        """
}

workflow {
    reads_ch = channel.fromPath("${params.input_dir}/*.fastq")

    trimmed_ch = trimming(reads_ch)
    fastqc_ch = fastqc_trimmed(trimmed_ch)
    aligned_ch = alignment(trimmed_ch)
    counts_ch = quantification(aligned_ch)

    all_results = counts_ch.mix(fastqc_ch)
    multiqc(all_results)
}
