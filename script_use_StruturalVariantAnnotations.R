##########################################################
##########################################################

# CRAN packages
library(devtools)
library(stringr)
library(dplyr)
library(ggplot2)

# bioconductor packages
library(VariantAnnotation)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

##########################################################
##########################################################
library(StructuralVariantAnnotation) 
# install_github("d-cameron/StructuralVariantAnnotation")
##########################################################
##########################################################

# a <- "./vcf.DELLY.SPCG-EWS008_2D.vcf"
  a <- "./vcf.DELLY.SPCG-EWS008_13R.vcf" 

vcf <- readVcf(a, "hg38")

rowRanges(vcf)
info(vcf)
info(header(vcf))
geno(vcf)
geno(header(vcf))

##########################################################
# breakpointRanges
# extracts the structural variants as a BreakendGRanges
##########################################################
# convert to breakend GRanges
##########################################################
##########################################################

gr <- breakpointRanges(vcf)

# retain only primary chromosomes
seqlevelsStyle(gr) <- "UCSC"
gr <- gr[seqnames(gr) %in% paste0("chr", c(1:22, "X", "Y"))]
seqlevels(gr) <- paste0("chr", c(1:22, "X", "Y"))
gr <- breakpointRanges(vcf)

##########################################################
##########################################################

gr
# remove breakends that now don't have a partner (eg: chr1 -> chrMT)
gr <- gr[gr$partner %in% names(gr)]

# annotate breakends with gene names and gene orientation

gns <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

hits <- as.data.frame(findOverlaps(gr, gns, ignore.strand=TRUE))

hits$SYMBOL <- biomaRt::select(org.Hs.eg.db, gns[hits$subjectHits]$gene_id, "SYMBOL")$SYMBOL

hits$gene_strand <- as.character(strand(gns[hits$subjectHits]))

hits <- hits %>%
  group_by(queryHits) %>%
  summarise(SYMBOL=paste(SYMBOL, collapse=","), gene_strand=paste0(gene_strand, collapse=""))

gr$SYMBOL <- ""
gr$geneStrand <- ""
gr$SYMBOL[hits$queryHits] <- hits$SYMBOL
gr$geneStrand[hits$queryHits] <- hits$gene_strand

# require the breakpoint to be between different genes
gr <- gr[gr$SYMBOL != partner(gr)$SYMBOL,]
gr <- gr[gr$SYMBOL != "" & partner(gr)$SYMBOL != "",]

##########################################################
##########################################################
##########################################################
# requires the breakpoint to possibly generate a fusion transcript
##########################################################
##########################################################
##########################################################

gr$couldBeThreePrimeStart <- str_detect(gr$geneStrand, stringr::fixed(as.character(strand(gr))))

gr$couldBeFivePrimeEnd <- str_detect(gr$geneStrand, stringr::fixed(ifelse(as.character(strand(gr))=="+", "-", "+")))

gr <- gr[(gr$couldBeThreePrimeStart & partner(gr)$couldBeFivePrimeEnd) |
         (gr$couldBeFivePrimeEnd & partner(gr)$couldBeThreePrimeStart),]

##########################################################
##########################################################
########################################################## writing the output :

write.table(gr, file="output_StructuralVariantAnnotation", 
                sep="\t", row.names=F)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


