params.index = 'samples.csv.orig'

Channel
    .fromPath(params.index)
    .splitCsv(header:true, quote:'\"')
    .map {row-> tuple(row.patient, row.samplename, file(row.bamfile), row.type) }
    .into { size_select_ch; readcounter_ch; ichorcna_ch}

process size_selection {

    cpus 8

    publishDir "${params.outdir}/${patient}/bams/", mode: "copy"
    input:
    set patient, sample, file(bamfile), type from size_select_ch

    output:
    tuple file("${sample}_${type}_ready.bam"), file("${sample}_${type}_ready.bam.bai"), sample, type, patient into processed_bam_ch

    script:
    dest = "${sample}_${type}_ready.bam"
    if (type == "cfdna") {
    """
    ${params.bam_size_select} single --index-source $bamfile ${dest} -c ${task.cpus}
    """
    } else {
    """
    cp -a $bamfile ${dest}
    samtools index -@${task.cpus} ${dest}
    """
    }

}

process read_counter {

    cpus 8

    publishDir "${params.outdir}/${patient}/wigs/", mode: "copy"
    input:
    set file(bamfile), file(bamindex), sample, type, patient from processed_bam_ch
    output:
    tuple file("${sample}.${type}.bin${params.readcounter_binsize}.wig"), sample, type, patient into wigfiles_ch
    script:

    """
    ${params.readcounter_path}  \
        --chromosome ${params.readcounter_chrs} \
        --window ${params.readcounter_binsize} \
        --quality ${params.readcounter_quality} \
        $bamfile \
        > ${sample}.${type}.bin${params.readcounter_binsize}.wig
    """

}

process ichorCNA {

    publishDir "${params.outdir}/${patient}/", mode: "copy"
    input:
    set file(wigfile), sample, type, patient from wigfiles_ch

    output:
    file "*.seg.txt"
    file "*.correctedDepth.txt"
    file "*.params.txt" into postprocess_ch
    file "*.seg"
    file "*.RData"
    path "${sample}"

    script:

    ploidy = "c(2,3)"
    compute_sc = "TRUE"
    compute_ploidy = "TRUE"
    sc_states = "c(1,3)"
    normal = "c(0.5,0.6,0.7,0.8,0.9,0.95)"

    if (type == "cfdna") {
        ploidy = "c(2)"
        compute_sc = "FALSE"
        compute_ploidy = "FALSE"
        sc_states = "c()"
        normal = "c(0.90, 0.98, 0.99, 0.995, 0.999)"
    }

    """
    mkdir -p ${sample}
    /usr/bin/Rscript ${params.ichorcna_script} \
            --id $sample \
            --WIG $wigfile \
            --gcWig ${params.gcwig} \
            --mapWig ${params.mapwig} \
            --centromere ${params.centromere} \
            --normalPanel ${params.normal_panel} \
            --estimatePloidy ${compute_ploidy} \
            --estimateScPrevalence ${compute_sc} \
            --estimateNormal ${params.estimate_normal} \
            --libdir "${params.ichorcna_libdir}" \
            --normal "${normal}" \
            --maxCN ${params.max_cn} \
            --scStates "${sc_states}" \
            --txnE ${params.txe} \
            --txnStrength ${params.tx_strength} \
            --minMapScore ${params.min_map_score} \
            --fracReadsInChrYForMale ${params.fraction_reads_male} \
            --maxFracGenomeSubclone ${params.max_frac_genome_subclone} \
            --maxFracCNASubclone ${params.max_frac_cna_subclone} \
            --includeHOMD ${params.include_homd} \
            --chrs "${params.chrs_to_use}" \
            --chrTrain "${params.chrs_to_train}" \
            --genomeStyle "${params.genome_style}" \
            --plotFileType "${params.plotfiletype}" \
            --plotYLim "${params.plotylim}" \
            --outDir ./
    """
}

process aggregate_table {

    publishDir "${params.outdir}/"
    input:
    file(all_params: "*.params.txt") from postprocess_ch.collect()
    output:
    file("ploidy_summary.txt")
    script:
    """
    awk 'FNR==1 && NR!=1 { while (/^Sample/) getline; }
    FNR < 3 {print}' ${all_params} > ploidy_summary.txt
    """
}
