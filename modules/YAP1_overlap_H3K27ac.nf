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
