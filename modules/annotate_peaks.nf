process ANNOTATE_PEAKS {
    conda 'bedtools'

    input:
    path H3K27ac_bed
    path YAP1_bed

    script:
    """
    #!/usr/bin/env Rscript
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(rtracklayer)
    library("org.Hs.eg.db")

    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

    # read in the peaks for YAP1
    YAP1<- import("YAP1_peaks.bed", format = "BED")
    YAP1

    YAP1_anno<- annotatePeak(YAP1, tssRegion=c(-3000, 3000),
                         TxDb=txdb, level = "gene", annoDb="org.Hs.eg.db",
                         sameStrand = FALSE, ignoreOverlap = FALSE,
                         overlap = "TSS")

    # some nice visualization you can do
    plotAnnoPie(YAP1_anno)
    upsetplot(YAP1_anno, vennpie=TRUE)

    # check the annotation
    head(as.data.frame(YAP1_anno))

    # you can save it to a txt file to your computer
    write.table(as.data.frame(YAP1_anno), "YAP1_peaks_anno.txt", row.names =F, col.names=T, sep ="\t", quote = F)
    """
}
