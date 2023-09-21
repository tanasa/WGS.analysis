
library("ggplot2")
library("reshape2")
library("limma")

a <- read.delim("file.DATAFRAME.samples.and.mutations.from.ANNOVAR.table.and.SUMMARY.SNV.and.CNV.txt", 
                 sep="\t", header=T, stringsAsFactors=T)

##############################################################
##############################################################

head(a)
dim(a)

##############################################################
##############################################################

plotWithHighlights(a$SNV.exonic, 
                   a$CNV.exonic,
                   status=a$Cancer_Subtype, 
                   bg.col="grey", 
                   # xlim=c(0,30), 
                   # ylim=c(-5,4), 
                   hl.cex=0.6, 
                   cex.main=0.8, 
                   cex.lab =0.8,
                   xlab="", ylab="", legend= "topright",
                   main="" )

##############################################################
##############################################################

ggplot(a, aes(x=SNV.exonic, y=CNV.exonic, col=Cancer_Subtype))+
       geom_point(size=2, aes(color=Cancer_Subtype)) +
       theme_bw() +
       scale_colour_manual(values = c( "HBL" = "#a6cee3",
                                      "RARE/OTHER" = "#1f78b4",
                                      "HCC" = "#b2df8a",
                                      "CNS" = "#33a02c",
                                      "WILMS" = "#fb9a99", 
                                      "NRSTS" = "#e31a1c",
                                      "RMS" = "#fdbf6f",
                                      "EWS" = "#ff7f00", 
                                      "OS" = "#cab2d6", 
                                      "LEUKEMIA/LYMPHOMA" = "#6a3d9a",
                                      "NBL" = "#ffff99" )) +
      ggtitle(expression(atop(bold("exonic CNV (cn.mops) vs exonic SNV+INDEL (MUTECT2)")))) +
      ylab("exonic CNV (cn.mops)") +
      xlab("exonic SNV (MUTECT2)") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(plot.title = element_text(size=12)) +
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank() 
            #, panel.border = element_blank() 
            ) 

##############################################################
##############################################################
ggsave("scatter.plot.EXONIC.SNV.vs.EXONIC.CNV.pdf", 
       height=20, width=20, units="cm")
##############################################################
##############################################################
