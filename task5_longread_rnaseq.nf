#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
 * =========================
 * Parameters
 * =========================
 */

params.reads = null
params.refFa = null
params.refGtf = null

params.outdir = "results"

/*
 * annotated / no_annotation
 */
params.scenario = "annotated"

/*
 * true  = annotation-guided Bambu
 * false = Bambu without annotations
 */
params.with_annotations = true



/*
 * =========================
 * 1. Minimap2 alignment
 * =========================
 */

process MINIMAP2_ALIGN {

    tag "${sample_id}"

    cpus 6
    memory '24 GB'

    input:
    tuple val(sample_id), path(reads)
    path refFa

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:

    /*
     * Direct RNA and cDNA require different
     * minimap2 parameters.
     */
    def mm2_opts = sample_id.contains("directRNA") \
        ? "-ax splice -uf -k14" \
        : "-ax splice"

    """
    minimap2 \
        -t ${task.cpus} \
        ${mm2_opts} \
        ${refFa} \
        ${reads} \
        > ${sample_id}.sam
    """
}



/*
 * =========================
 * 2. SAM -> sorted BAM
 * =========================
 */

process SAM_TO_BAM {

    tag "${sample_id}"

    cpus 4
    memory '8 GB'

    publishDir "${params.outdir}/bam",
        mode: 'copy'

    input:
    tuple val(sample_id), path(reads_sam)

    output:
    tuple val(sample_id),
          path("${sample_id}.bam"),
          path("${sample_id}.bam.bai")

    script:
    """
    samtools sort \
        -@ ${task.cpus} \
        ${reads_sam} \
        -o ${sample_id}.bam

    samtools index \
        -@ ${task.cpus} \
        ${sample_id}.bam
    """
}



/*
 * =========================
 * 3. QC
 * =========================
 */

process QC_BAM {

    tag "${sample_id}"

    cpus 2
    memory '4 GB'

    publishDir "${params.outdir}/qc",
        mode: 'copy'

    input:
    tuple val(sample_id),
          path(reads_bam),
          path(reads_bai)

    output:
    path("${sample_id}_qc.tsv")
    path("${sample_id}_flagstat.txt")
    path("${sample_id}_stats.txt")

    script:
    """
    samtools flagstat ${reads_bam} \
        > ${sample_id}_flagstat.txt

    samtools stats ${reads_bam} \
        > ${sample_id}_stats.txt


    # Primary reads:
    # exclude secondary (256) and supplementary (2048)
    primary_total=\$(samtools view -c -F 2304 ${reads_bam})

    # Primary mapped reads:
    # additionally exclude unmapped (4)
    primary_mapped=\$(samtools view -c -F 2308 ${reads_bam})


    mapping_pct=\$(awk \
        -v m=\$primary_mapped \
        -v t=\$primary_total \
        'BEGIN {
            if (t > 0)
                printf "%.2f", 100*m/t;
            else
                printf "0.00"
        }')


    avg_length=\$(awk -F '\\t' \
        '\$1=="SN" && \$2=="average length:" {print \$3}' \
        ${sample_id}_stats.txt)


    # User-defined QC criteria:
    # mapping rate >= 70%
    # at least 1000 primary reads

    pass_mapping=\$(awk \
        -v p=\$mapping_pct \
        'BEGIN {print (p >= 70) ? 1 : 0}')


    if [ "\$primary_total" -ge 1000 ] \
       && [ "\$pass_mapping" -eq 1 ]; then

        qc_status="PASS"

    else

        qc_status="FAIL"

    fi


    printf "sample\\tprimary_reads\\tprimary_mapped\\tmapping_pct\\taverage_length\\tQC\\n" \
        > ${sample_id}_qc.tsv

    printf "${sample_id}\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \
        "\$primary_total" \
        "\$primary_mapped" \
        "\$mapping_pct" \
        "\$avg_length" \
        "\$qc_status" \
        >> ${sample_id}_qc.tsv
    """
}



/*
 * =========================
 * 4. Bambu
 * =========================
 */

process BAMBU {

    tag "${params.scenario}"

    cpus 6
    memory '48 GB'

    publishDir "${params.outdir}/${params.scenario}/bambu",
        mode: 'copy'

    input:
    path refFa
    path refGtf
    path bam_bundle
    val with_annotations

    output:
    path "counts_transcript.txt"
    path "counts_gene.txt"
    path "extended_annotations.gtf"
    path "bambu_se.rds"

    script:

    /*
     * Scenario 1
     */
    def bambu_call_with_annotation = """
    annotations <- prepareAnnotations("${refGtf}")

    se <- bambu(
        reads = bam_files,
        annotations = annotations,
        genome = "${refFa}",
        ncore = ${task.cpus}
    )
    """


    /*
     * Scenario 2
     * Follows the workshop approach for
     * transcript discovery without annotation.
     */
    def bambu_call_without_annotation = """
    se <- bambu(
        reads = bam_files,
        genome = "${refFa}",
        NDR = 1,
        opt.discovery = list(
            min.readFractionByEqClass = 0.2
        ),
        ncore = ${task.cpus}
    )
    """


    def bambu_call = with_annotations \
        ? bambu_call_with_annotation \
        : bambu_call_without_annotation


    """
    #!/usr/bin/env Rscript --vanilla

    library(bambu)

    bam_files <- Sys.glob("*.bam")

    ${bambu_call}

    writeBambuOutput(
        se,
        path = "./"
    )

    saveRDS(
        se,
        file = "bambu_se.rds"
    )
    """
}



/*
 * =========================
 * Workflow
 * =========================
 */

workflow {

    /*
     * FASTQ channel
     */

    reads_ch = Channel
        .fromPath(params.reads, checkIfExists: true)
        .map { fq ->

            def sid = fq.name
                .replaceFirst(/\\.fastq\\.gz$/, '')
                .replaceFirst(/\\.fq\\.gz$/, '')

            tuple(sid, fq)
        }


    /*
     * Reference files are value channels,
     * so they can be reused for every sample.
     */

    refFa_ch = Channel.value(
        file(params.refFa)
    )

    refGtf_ch = Channel.value(
        file(params.refGtf)
    )


    /*
     * Alignment
     */

    MINIMAP2_ALIGN(
        reads_ch,
        refFa_ch
    )


    /*
     * BAM conversion
     */

    SAM_TO_BAM(
        MINIMAP2_ALIGN.out
    )


    /*
     * QC
     */

    QC_BAM(
        SAM_TO_BAM.out
    )


    /*
     * Collect BAM + BAI files so Bambu receives
     * all samples together.
     */

    bam_bundle_ch = SAM_TO_BAM.out
        .flatMap { sid, bam, bai ->
            [bam, bai]
        }
        .collect()


    annotation_mode_ch = Channel.value(
        params.with_annotations
            .toString()
            .toBoolean()
    )


    /*
     * Bambu
     */

    BAMBU(
        refFa_ch,
        refGtf_ch,
        bam_bundle_ch,
        annotation_mode_ch
    )
}
