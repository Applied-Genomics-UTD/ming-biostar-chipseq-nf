process ANNOTATE_PEAKS {
    conda './envs/annotate_peaks.yml'

    input:
    path H3K27ac_bed
    path YAP1_bed

    output:
    path "YAP1_peaks_anno.txt"
    path "Rplots.pdf"
    path "YAP_KEGG_pathway_genes.txt"

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

    # https://www.biostarhandbook.com/chip-seq-downstream-analysis-1.html#how-do-i-do-pathway-enrichment-analysis-for-the-peaks
    library(clusterProfiler)

    ## GO term enrichment
    ego <- enrichGO(gene           = as.data.frame(YAP1_anno)$SYMBOL,
                    OrgDb         = org.Hs.eg.db,
                    keytype       = "SYMBOL",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
    # visualization
    dotplot(ego, showCategory = 20)


    ## Kegg pathway enrichment, need the Entrez ID
    kk<- enrichKEGG(gene    = as.data.frame(YAP1_anno)$geneId,
            organism     = 'hsa',
            pvalueCutoff = 0.05)

    dotplot(kk, showCategory = 20, title = "YAP1 binding pathway enrichment")

    ## you can write the result to a tsv file
    write.table(kk@result, "YAP_KEGG_pathway_genes.txt", sep = "\t", col.names = T, row.names = F, quote =F)
    """
}
