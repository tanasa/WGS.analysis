############################################################
############################################################
############################################################
############################################################

suppressMessages(library(ShatterSeek))
suppressMessages(library(dplyr)) 
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(tidyverse))

############################################################
############################################################
############################################################
############################################################
# to read the DATAFRAME with SV
# to read the DATAFRAME with CNV
# to read the NAME of the SAMPLE
############################################################
############################################################
############################################################
############################################################ to use the R script
# we use the R script in the following way :

# sv="${line}.manta.somatic.svs.hg38.v2.bedpe"
# cna="${line}_allelic_states.txt"
# sample="$line"

# echo $sv
# echo $cna
# echo $sample

# Rscript --vanilla  \
# the_script_to_compute_chromothripsis.R \
# $sv \
# $cna \
# $sample

############################################################
############################################################
############################################################
############################################################ SV file is in BEDPE format
# the FORMAT of the SV file is : 
# 1  866187   866197   7  96949165 96949165 MantaBND:0:821092:821096:0:0:0:0 ,  +  +  TRA   1  866187   7  96949166
# 1  1329165  1329174  18 45390773 45390773 MantaBND:0:1040335:1268996:0:5:0:1  ,  +  +  TRA   1  1329165  18 45390774
# 1  1613901  1613901  4  163533566   163533566   MantaBND:0:186358:186369:0:0:0:1 ,  +  +  TRA   1  1613902  4  163533567
############################################################
############################################################
############################################################
############################################################ CNV file contains the ALLELIC states
# the FORMAT of the CNV file is : 
# Sample   Chromosome  Start End   Copy_Nr_Raw CopyNr   A  B  LOH   Theta Theta_Exp   n_SNPs   Is_Subclonal_CN   Subclonal_P_value Is_Inconsistent_State
# EA5040686   chr1  10001 8881929  2.07513  2  1  1  0  0.275862 0.254501 7787  0  0.15442  0
# EA5040686   chr1  8881930  10594318 2.16317  2  1  1  0  0.318182 0.254501 1487  1  2.18195e-05 1
# EA5040686   chr1  10594319 15950516 2.0449   2  1  1  0  0.272727 0.254501 5511  0  0.224326 0
############################################################
############################################################
############################################################
############################################################

args = commandArgs(trailingOnly=TRUE)

SV= args[1]
CNA = args[2]
SAMPLE = args[3]

############################################################
############################################################3#

SAMPLE = paste(SAMPLE, "shatterseek", sep=".")
sv.df = read.table(SV, header=F, sep="\t", stringsAsFactors=F)
cnv.df = read.table(CNA, header=T, sep="\t", stringsAsFactors=F)

############################################################
############################################################
# colnames(cnv.df)
# [1] "Sample"                "Chromosome"            "Start"
# [4] "End"                   "Copy_Nr_Raw"           "CopyNr"
# [7] "A"                     "B"                     "LOH"
# [10] "Theta"                 "Theta_Exp"             "n_SNPs"
# [13] "Is_Subclonal_CN"       "Subclonal_P_value"     "Is_Inconsistent_State"
############################################################
############################################################

# mainChr = c(as.character(1:22),'x','X')
mainChr = c(as.character(1:22))
sv.df = sv.df[sv.df$V1 %in% mainChr,]
sv.df = sv.df[sv.df$V4 %in% mainChr,]

# mainChr = c(as.character(1:22),'x','X')
mainChr = c(as.character(1:22))
mainChr.cnv = paste("chr",mainChr, sep="")
cnv.df = dplyr::filter(cnv.df, Chromosome %in% mainChr.cnv)

# sv.df$V1 = paste("chr", sv.df$V1, sep="")
# sv.df$V4 = paste("chr", sv.df$V4, sep="")

SV_data <- SVs(chrom1=as.character(sv.df$V1), 
			   pos1=as.numeric(sv.df$V2),
			   chrom2=as.character(sv.df$V4), 
			   pos2=as.numeric(sv.df$V5),
			   SVtype=as.character(sv.df$V11), 
			   strand1=as.character(sv.df$V9),
			   strand2=as.character(sv.df$V10))


# to change the names of the chromosomes
cnv.df$Chr = str_replace(cnv.df$Chromosome, "chr","")

CN_data <- CNVsegs(chrom=as.character(cnv.df$Chr),
				   start=cnv.df$Start,
				   end=cnv.df$End,
				   total_cn=cnv.df$CopyNr)

chromothripsis <- shatterseek(
                SV.sample=SV_data,
                seg.sample=CN_data,
                genome="hg38")              

############################################################
############################################################

# str(SV_data)
# str(CN_data)
# str(chromothripsis)

chromosomes = unique(cnv.df$Chr)
length(chromosomes)
head(chromothripsis@chromSummary, 2)

columns_of_interest <- "start"
chromo_summary.df = as.data.frame(chromothripsis@chromSummary)
v
head(chromo_summary.df, 2)
dim(chromo_summary.df)

chromo_summary.df.filtered <- chromo_summary.df[complete.cases(chromo_summary.df[, columns_of_interest]), ]

head(chromo_summary.df.filtered, 2)
dim(chromo_summary.df.filtered)

chromosomes.filtered = chromo_summary.df.filtered$chrom
str(chromosomes.filtered)
chromosomes.filtered
length(chromosomes.filtered)

############################################################
############################################################

for (i in 1:length(chromosomes.filtered) ){

try(
{

chr = chromosomes.filtered[[i]]

plots_chr <- plot_chromothripsis(ShatterSeek_output = chromothripsis, 
               chr = chr, 
               sample_name=SAMPLE, genome="hg38")

plots_chr = arrangeGrob(plots_chr[[1]],
                         plots_chr[[2]],
                         plots_chr[[3]],
                         plots_chr[[4]],
                         nrow=4,ncol=1,heights=c(0.2,.4,.4,.4))  


# pdf(file = paste(SAMPLE,"chromothripsis", chr, "pdf", sep="."),
#         width = 8, height = 10, pointsize = 20)
# plot_grid(plots_chr)
# dev.off()

# png(filename = paste(SAMPLE,"chromothripsis", chr, "png", sep="."),
#         width = 600, height = 800, units = "px", pointsize = 10)
# plot_grid(plots_chr)
# dev.off()

 ggsave(plot_grid(plots_chr),
    filename = paste(SAMPLE,"chromothripsis", chr, "jpeg", sep="."),
    width = 40,
    height = 30,
    units = "cm"
 )

}, silent = FALSE)

}

# chromothripsis@detail
# chromothripsis@detail$SV
# chromothripsis@detail$graph
# chromothripsis$SVinter
# chromothripsis$CNV
# str(chromothripsis)

# chromo_summary = as.data.frame(chromothripsis@chromSummary)

############################################################
############################################################

# to print the results

write.table(chromo_summary.df, 
            file=paste(SAMPLE,"chromothripsis.SUMMARY.txt", sep="."),
            quote=FALSE, sep="\t", col.names = TRUE, row.names = FALSE)

write.table(chromo_summary.df.filtered, 
            file=paste(SAMPLE,"chromothripsis.SUMMARY.filtered.chromo.txt", sep="."),
            quote=FALSE, sep="\t", col.names = TRUE, row.names = FALSE)
            
            
write.table(chromothripsis@detail$num.chromth, 
            file=paste(SAMPLE,"chromothripsis.SUMMARY.number.chromo.txt", sep="."),
            quote=FALSE, sep="\t", col.names = TRUE, row.names = FALSE)

write.table(chromothripsis@detail$maxSVs, 
            file=paste(SAMPLE,"chromothripsis.SUMMARY.max.SV.txt", sep="."),
            quote=FALSE, sep="\t", col.names = TRUE, row.names = FALSE)

cluster_size = as.data.frame(chromothripsis@detail$maxClusterSize)

write.table(cluster_size, 
            file=paste(SAMPLE,"chromothripsis.SUMMARY.max.Cluster.Size.txt", sep="."),
            quote=FALSE, sep="\t", col.names = TRUE, row.names = FALSE)
         
############################################################
############################################################ tests
 
# i = 22
# chr = chromosomes[[i]]

# plots_chr <- plot_chromothripsis(ShatterSeek_output = chromothripsis, 
#               chr = chr, 
#               sample_name=SAMPLE, genome="hg38")

# plots_chr = arrangeGrob(plots_chr[[1]],
#                         plots_chr[[2]],
#                         plots_chr[[3]],
#                         plots_chr[[4]],
#                         nrow=4,ncol=1,heights=c(0.2,.4,.4,.4))  


### pdf(file = paste(SAMPLE,"chromothripsis", chr, "pdf", sep="."),
###         width = 8, height = 10, pointsize = 20)
### plot_grid(plots_chr)
### dev.off()

# png(filename = paste(SAMPLE,"chromothripsis", chr, "png", sep="."),
#         width = 600, height = 800, units = "px", pointsize = 10)
# plot_grid(plots_chr)
# dev.off()

############################################################
############################################################
