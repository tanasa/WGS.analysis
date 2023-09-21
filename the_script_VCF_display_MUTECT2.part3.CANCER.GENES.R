library("ggplot2")

############################################################

x <- read.delim("SUMMARY_of_SAMPLES.19july2017.regarding.number.INDELS.and.SNV..and.from.ANNOVAR.table.protein.coding.regions.and.general.and.in.CANCER.GENES.info.for.R",
header=T, stringsAsFactors=F, sep="\t")

head(x)

########################################################
########################################################

# Turn your 'INTERESTING' column into a character vector
x$SAMPLE <- as.character(x$SAMPLE)
# Then turn it back into an ordered factor
x$SAMPLE <- factor(x$SAMPLE, levels=unique(x$SAMPLE))

## doing the same for the field "Cancer_Subtype"

# Turn your 'INTERESTING' column into a character vector
x$Cancer_Subtype <- as.character(x$Cancer_Subtype)
# Then turn it back into an ordered factor
x$Cancer_Subtype <- factor(x$Cancer_Subtype, levels=unique(x$Cancer_Subtype))

## here to compute the MUTATIONS that are EXONIC (exonic, ncRNA exonic), and alter SPLICING : 
## only nonSYN SNV and INDELS
## in the CANCER GENES

########################################################
########################################################

x$cancer_genes <- x$frameshift_deletion  + 
                  x$frameshift_insertion  +
                  x$nonframeshift_deletion  +
                  x$nonframeshift_insertion  +
                  x$nonsyn_SNV  +
                  x$stopgain + 
                  x$other_mutations

head(x)

######################################################
######################################################

# pdf("MUTECT2 predicted mutations.pdf")

ggplot(x, aes(y=cancer_genes, x=Cancer_Subtype)) +
      # geom_boxplot() + 
      geom_point(size=4, aes(shape=Treatment_Status, color=Cancer_Subtype)) + 
      # geom_point(size=4, aes(color=Cancer_Subtype)) + 
      scale_shape_manual(values=c(15,16,17,18)) +
      theme_bw() +
      ggtitle("mutations (nonsyn-SNV, INDEL, splicing) in CANCER GENES per CANCER SUBTYPE") +
      ylab("mutations (nonsyn-SNV, INDEL, splicing) in CANCER GENES") +
      xlab("cancer subtype") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(plot.title = element_text(size=12)) +
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            #panel.border = element_blank()
            )

## theme(text=element_text(family="Arial", size=14))
## dev.off()

ggsave("figure -- MUTECT2 mutations (nonsyn-SNV, INDEL, and SPLICING) in CANCER GENES per cancer subtype.pdf")

###########################################################
###########################################################
###########################################################
