#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.blacklist_bed = file("https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz")
params.H3K27ac_peaks = file("/scratch/applied-genomics/chipseq/ming-results/bwa/mergedLibrary/macs2/broadPeak/WT_H3K27ac_peaks.broadPeak")
params.YAP1_peaks = file("/scratch/applied-genomics/chipseq/ming-results/bwa/mergedLibrary/macs2/broadPeak/WT_YAP1_peaks.broadPeak")

process compare_peak_sets {
    conda 'bedtools'

    input:
    path H3K27ac
    path YAP1
    path blacklist

    output:
    path 'H3K27ac_filtered_peaks.bed', emit: H3K27ac
    path 'YAP1_filtered_peaks.bed', emit: YAP1

    script:
    """
    # How many H3K27ac peaks overlap with black-listed regions?
    bedtools intersect -a $H3K27ac -b $blacklist -wa | wc -l
    #14

    # exclude those peaks
    bedtools intersect -a $H3K27ac -b $blacklist -v > H3K27ac_filtered_peaks.bed

    # do the same for YAP1
    bedtools intersect -a $YAP1 -b $blacklist -v > YAP1_filtered_peaks.bed
    """
}

process YAP1_overlap_H3K27ac {
    conda 'bedtools'

    input:
    path H3K27ac_bed
    path YAP1_bed

    script:
    """
    bedtools intersect -a $YAP1_bed -b $H3K27ac_bed -wa | wc -l
    #1882

    bedtools intersect -a $YAP1_bed -b $H3K27ac_bed -wa | sort | uniq | wc -l
    #1882

    bedtools intersect -a $H3K27ac_bed -b $YAP1_bed -wa | wc -l
    #1882

    bedtools intersect -a $H3K27ac_bed -b $YAP1_bed -wa | sort | uniq | wc -l
    #1772

    bedtools intersect -a H3K27ac_filtered_peaks.bed -b YAP1_filtered_peaks.bed -wa | sort | uniq -c | sort -k1,1nr | head
    """
}


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

}
