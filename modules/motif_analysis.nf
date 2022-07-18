process MOTIF_ANALYSIS {
    conda 'bedtools'
    publishDir 'results/motif_analysis'

    input:
    path H3K27ac_bed
    path YAP1_summits

    output:
    path "UCSC_hg19_genome.fa.gz"
    path "YAP1_500bp.fa"

    script:
    """
    cat $YAP1_summits | awk '\$2=\$2-249, \$3=\$3+250' OFS="\t" > YAP1_500bp_summits.bed

    cat *fa.gz > UCSC_hg19_genome.fa.gz
    gunzip UCSC_hg19_genome.fa.gz

    # Use betools get fasta http://bedtools.readthedocs.org/en/latest/content/tools/getfasta.html
    bedtools getfasta -fi UCSC_hg19_genome.fa -bed YAP1_500bp_summits.bed -fo YAP1_500bp.fa
    """

}
