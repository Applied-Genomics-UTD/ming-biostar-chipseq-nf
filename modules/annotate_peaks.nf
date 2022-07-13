process ANNOTATE_PEAKS {
    conda './envs/annotate_peaks.yml'

    input:
    path H3K27ac_bed
    path YAP1_bed

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(ChIPseeker))
    suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
    suppressPackageStartupMessages(library(rtracklayer))
    suppressPackageStartupMessages(library("org.Hs.eg.db"))

    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

    # read in the peaks for YAP1
    YAP1 <- readPeakFile("$YAP1_bed", as="GRanges")
    YAP1

    YAP1_anno<- annotatePeak(YAP1, tssRegion=c(-3000, 3000),
                         TxDb=txdb, level = "gene", annoDb="org.Hs.eg.db",
                         sameStrand = FALSE, ignoreOverlap = FALSE,
                         overlap = "TSS")

    # some nice visualization you can do
    plotAnnoPie(YAP1_anno)
    upsetplot(YAP1_anno, vennpie=FALSE)

    # check the annotation
    head(as.data.frame(YAP1_anno))

    # you can save it to a txt file to your computer
    write.table(as.data.frame(YAP1_anno), "YAP1_peaks_anno.txt", row.names =F, col.names=T, sep ="\t", quote = F)
    """
}
