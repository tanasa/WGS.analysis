library("ggplot2")

############################################################

x <- read.delim("SUMMARY_of_SAMPLES.19july2017.regarding.number.INDELS.and.SNV.and.from.ANNOVAR.table.general.distribution.INDEL.SNV.info.for.R",
header=T, stringsAsFactors=F, sep="\t")

head(x)

########################################################
########################################################

# Turn your 'INTERESTING' column into a character vector
x$SAMPLE <- as.character(x$SAMPLE)
# Then turn it back into an ordered factor
x$SAMPLE <- factor(x$SAMPLE, levels=unique(x$SAMPLE))

# doing the same for the field "Cancer_Subtype"

# Turn your 'INTERESTING' column into a character vector
x$Cancer_Subtype <- as.character(x$Cancer_Subtype)
# Then turn it back into an ordered factor
x$Cancer_Subtype <- factor(x$Cancer_Subtype, levels=unique(x$Cancer_Subtype))

## to add "shape=TUMOR_STAGE" in the NEXT VERSION

######################################################
######################################################

## pdf("figure -- MUTECT2 (filtered) mutations (SNV+INDEL) per cancer subtype.pdf", height=10, width=10)

ggplot(x, aes(y=TOTAL_f, x=Cancer_Subtype)) +
      # geom_boxplot() + 
      geom_point(size=4, aes(shape=Treatment_Status, color=Cancer_Subtype)) + 
      # geom_point(size=4, aes(color=Cancer_Subtype)) + 
      scale_shape_manual(values=c(15,16,17,18)) +
      theme_bw() +
      ggtitle("MUTATIONS (MUTECT2) per CANCER SUBTYPE") +
      ylab("mutations (SNV + INDEL)") +
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

ggsave("figure -- MUTECT2 (filtered) mutations (SNV+INDEL) per cancer subtype.pdf")

###########################################################
###########################################################
###########################################################

## pdf("MUTECT2 predicted SNV.pdf")

ggplot(x, aes(y=SNV_f, x=Cancer_Subtype)) +
      ## geom_boxplot() + 
      geom_point(size=4, aes(shape=Treatment_Status, color=Cancer_Subtype)) + 
      # geom_point(size=4, aes(color=Cancer_Subtype)) + 
      scale_shape_manual(values=c(15,16,17,18)) +
      theme_bw() +
      ggtitle("SNV (MUTECT2) per CANCER SUBTYPE") +
      ylab("SNV") +
      xlab("cancer subtype") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(plot.title = element_text(size=12)) +
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            #panel.border = element_blank()
            )


ggsave("figure -- MUTECT2 SNV per cancer subtype.pdf")

######################################################
#####################################################

## pdf("MUTECT2 predicted insertions.pdf")

ggplot(x, aes(y=INSERTION_f, x=Cancer_Subtype)) +
      # geom_boxplot() + 
      geom_point(size=4, aes(shape=Treatment_Status, color=Cancer_Subtype)) + 
      # geom_point(size=4, aes(color=Cancer_Subtype)) + 
      scale_shape_manual(values=c(15,16,17,18)) +
      theme_bw() +
      ggtitle("INSERTIONS (MUTECT2) per CANCER SUBTYPE") +
      ylab("INS") +
      xlab("cancer subtype") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(plot.title = element_text(size=12)) +
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            #panel.border = element_blank()
            )

ggsave("figure -- MUTECT2 INSERTIONS per cancer subtype.pdf")

## pdf("MUTECT2 predicted deletions.pdf")

ggplot(x, aes(y=DELETION_f, x=Cancer_Subtype)) +
      # geom_boxplot() + 
      geom_point(size=4, aes(shape=Treatment_Status, color=Cancer_Subtype)) + 
      # geom_point(size=4, aes(color=Cancer_Subtype)) + 
      scale_shape_manual(values=c(15,16,17,18)) +
      theme_bw() +
      ggtitle("DELETIONS (MUTECT2) per CANCER SUBTYPE") +
      ylab("DEL") +
      xlab("cancer subtype") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(plot.title = element_text(size=12)) +
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
            #panel.border = element_blank()
            )

## dev.off()
ggsave("figure -- MUTECT2 DELETIONS per cancer subtype.pdf")

#########################################################################
#########################################################################

# pdf("MUTECT2 predicted non-syn SNVs and INDELs.pdf")

# ggplot(x, aes(y=nonsynonymous_SNV_and_INDELS, x=TUMOR_TYPE)) +
#      geom_boxplot() + 
#      geom_point(size=4, aes(color=TUMOR_TYPE, shape=TUMOR_STAGE)) + 
#      scale_shape_manual(values=c(15,16,17,18)) +
#      theme_bw() +
#      ggtitle("MUTECT2 predicted nonsynonymous SNV and INDELS") +
#      ylab("nonsynonymous SNV and INDELS") +
#      xlab("cancer type") +
#      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#      theme(plot.title = element_text(size=12)) +
#      theme(plot.background = element_blank(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank()
#            #panel.border = element_blank()
#            )

# dev.off()
# ggsave("MUTECT2 predicted nonsyn SNVS and INDELS per cancer type v3.pdf")


###########################################################
###########################################################
###########################################################

### IN ORDER to remove the GRID from the PICTURES :

# + theme(axis.line = element_line(colour = "black"),
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(),
#    panel.border = element_blank(),
#    panel.background = element_blank()
#        )

###########################################################
###########################################################
###########################################################