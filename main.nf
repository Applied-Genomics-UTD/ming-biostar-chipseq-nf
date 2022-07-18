#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.blacklist_bed = file("https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz")
params.H3K27ac_peaks = file("/scratch/applied-genomics/chipseq/ming-results/bwa/mergedLibrary/macs2/broadPeak/WT_H3K27ac_peaks.broadPeak")
params.YAP1_peaks = file("/scratch/applied-genomics/chipseq/ming-results/bwa/mergedLibrary/macs2/broadPeak/WT_YAP1_peaks.broadPeak")
params.YAP1_summits = file("/scratch/applied-genomics/chipseq/ming-results/bwa/mergedLibrary/macs2/broadPeak/WT_YAP1_summits.bed")
// http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*
params.hg19_chrom = file("ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*")

include { compare_peak_sets } from "./modules/compare_peak_sets"
include { YAP1_overlap_H3K27ac } from "./modules/YAP1_overlap_H3K27ac.nf"
include { ANNOTATE_PEAKS } from "./modules/annotate_peaks.nf"
include { MOTIF_ANALYSIS } from "./modules/motif_analysis.nf"

workflow {
    compare_peak_sets (
        params.H3K27ac_peaks,
        params.YAP1_peaks,
        params.blacklist_bed
    )

    YAP1_overlap_H3K27ac (
        compare_peak_sets.out.H3K27ac,
        compare_peak_sets.out.YAP1
    )

    ANNOTATE_PEAKS (
        params.H3K27ac_peaks,
        params.YAP1_peaks
    )
    MOTIF_ANALYSIS (
        params.YAP1_summits,
        params.hg19_chrom
    )
}
