##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################

suppressMessages(require(sigrap))
suppressMessages(require(gpgr))
suppressMessages(require(MutationalPatterns))
suppressMessages(require(devtools))
suppressMessages(require(dplyr))
suppressMessages(require(patchwork))
suppressMessages(library(tidyverse))
ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
suppressMessages(require(ref_genome, character.only = TRUE))

##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################

# VCF = "SRR8652105_WGS-MCF7_BREAST.somatic.strelka2.only.mutations.hg38.vcf"
# VCF.NAME = basename(VCF)
# VCF.NAME <- gsub('.vcf', '', VCF.NAME)
# VCF.NAME

# Check if a command line argument is provided
if (length(commandArgs(trailingOnly = TRUE)) > 0) {

  VCF <- commandArgs(trailingOnly = TRUE)[1]
  VCF.NAME = basename(VCF)
  VCF.NAME <- gsub('.vcf', '', VCF.NAME)
  VCF.NAME
  
} else {
	
  cat("Usage: Rscript read_file.R <file_name>\n")
  stop("Script is stopping.")
  
}

##################################################################################################
##################################################################################################

tryCatch({
	
params <- list(
  vcf = VCF,
  nm = basename(VCF),
  outdir = paste(basename(VCF), "MutationalPatterns.outputs", sep="")
)

gr <- MutationalPatterns::read_vcfs_as_granges(
  vcf_files = params$vcf,
  sample_names = params$nm,
  genome = ref_genome,
  group = "auto+sex",
  type = "all"
)

head(gr, 2)
tail(gr, 2)

}, error = function(e)-999) 

################################################################################################### SNVS (SBS)
################################################################################################### PLOTS
###################################################################################################
################################################################################################### SNVS (SBS)

tryCatch({

snv_counts <- sigrap::sig_count_snv(vcf_gr = gr, ref_genome = ref_genome)

p_snvs <- sigrap::sig_plot_snv(
  gr_snv = snv_counts$gr_snv, snv_counts = snv_counts$snv_counts,
  ref_genome = ref_genome
)

# png(filename = paste(VCF.NAME ,"riverPlot.png", sep="."),
#     width = 1600, height = 600, units = "px", pointsize = 12)
# p_snvs$p_river + p_snvs$p_heatmap + p_snvs$p_spectrum + p_snvs$p_96_profile + patchwork::plot_layout(ncol = 1)
# dev.off()

p_snvs$p_river + p_snvs$p_heatmap + p_snvs$p_spectrum + p_snvs$p_96_profile + patchwork::plot_layout(ncol = 1)  
ggsave(filename = paste(VCF.NAME ,"SBS.riverPlot.png", sep="."), 
        width = 3400, height = 3400, units = "px", pointsize = 12)
    
    
}, error = function(e)-999)    

##################################################################################################
##################################################################################################    
##################################################################################################
##################################################################################################

# options(repr.plot.width = 16, repr.plot.height = 20)
# p_snvs$p_river + p_snvs$p_heatmap + p_snvs$p_spectrum + p_snvs$p_96_profile + patchwork::plot_layout(ncol = 1)

# options(repr.plot.width = 32, repr.plot.height = 8)
# p_snvs$p_river  + patchwork::plot_layout(ncol = 1)

# options(repr.plot.width = 32, repr.plot.height = 8)
# p_snvs$p_heatmap + patchwork::plot_layout(ncol = 1)

# options(repr.plot.width = 4, repr.plot.height = 6)
# p_snvs$p_spectrum + patchwork::plot_layout(ncol = 1)

# options(repr.plot.width = 22, repr.plot.height = 10)
# p_snvs$p_96_profile + patchwork::plot_layout(ncol = 1)

################################################################################################### SNVS (SBS)
################################################################################################### Signature Contributions
###################################################################################################
################################################################################################### SNVS (SBS)

tryCatch({

sigs_snv_2015 <-
  sigrap::cosmic_signatures_2015 |>
  {
    \(sigs) sigrap::sig_contribution(mut_mat = snv_counts$snv_counts, signatures = sigs)
  }()

    
sigs_snv_2015.df = as.data.frame(sigs_snv_2015)
sigs_snv_2015.df$sample = VCF.NAME
    
write.table(sigs_snv_2015.df, file=paste(VCF.NAME, "contributions.SIG.2015.txt", sep="."), 
sep="\t", col.names=TRUE,  row.names = FALSE, append = FALSE, quote = FALSE)

}, error = function(e)-999)
         
##################################################################################################
##################################################################################################
##################################################################################################
################################################################################################## 

tryCatch({
    
head(sigs_snv_2015.df, 2)
    
}, error = function(e)-999)

##################################################################################################
################################################################################################## 

tryCatch({
    
sigs_snv_2020 <-
  MutationalPatterns::get_known_signatures(
    muttype = "snv",
    incl_poss_artifacts = TRUE
  ) |>
  {
    \(sigs) sigrap::sig_contribution(mut_mat = snv_counts$snv_counts, signatures = sigs)
  }()
    
sigs_snv_2020.df = as.data.frame(sigs_snv_2020)
sigs_snv_2020.df$sample = VCF.NAME
    
write.table(sigs_snv_2020.df, file=paste(VCF.NAME, "contributions.SBS.2020.txt", sep="."), 
sep="\t", col.names=TRUE,  row.names = FALSE, append = FALSE, quote = FALSE)
    
}, error = function(e)-999)
         
##################################################################################################
##################################################################################################
##################################################################################################
################################################################################################## 

tryCatch({
    
head(sigs_snv_2020.df, 2)   
    
}, error = function(e)-999)

##################################################################################################
################################################################################################## 
################################################################################################## 2015
################################################################################################## HTML

tryCatch({
    
sigs_snv_2015 |>
  sigrap::sig_contribution_table(type = "Sig", outdir = params$outdir) |>
  knitr::kable(format = "html") |>
  kableExtra::kable_styling(c("hover", "striped"), font_size = 12) |>
  kableExtra::scroll_box(height = "400px")
    
}, error = function(e)-999)
         
##################################################################################################
##################################################################################################
################################################################################################## 2020
################################################################################################## HTML

tryCatch({
    
sigs_snv_2020 |>
  sigrap::sig_contribution_table(type = "SBS", outdir = params$outdir) |>
  knitr::kable(format = "html") |>
  kableExtra::kable_styling(c("hover", "striped"), font_size = 12) |>
  kableExtra::scroll_box(height = "400px")
    
}, error = function(e)-999)
         
###################################################################################################
###################################################################################################
################################################################################################### DBS
###################################################################################################
################################################################################################### PLOTS
################################################################################################### DBS 

tryCatch({
    
dbs_counts <- sigrap::sig_count_dbs(vcf_gr = gr)
p_dbs <- sigrap::sig_plot_dbs(dbs_counts = dbs_counts)

p_dbs$p_dbs_main / p_dbs$p_dbs_cont

ggsave(filename = paste(VCF.NAME ,"DBS.png", sep="."), 
        width = 3200, height = 2200, units = "px", pointsize = 12) 

}, error = function(e)-999)
         
##################################################################################################
##################################################################################################

# options(repr.plot.width = 12, repr.plot.height = 6)
# p_dbs$p_dbs_main 

# options(repr.plot.width = 22, repr.plot.height = 6)
# p_dbs$p_dbs_cont

# options(repr.plot.width = 22, repr.plot.height = 14)
# p_dbs$p_dbs_main / p_dbs$p_dbs_cont

##################################################################################################
################################################################################################## 

# tryCatch({
#
# p_dbs$p_dbs_main 
# ggsave(filename = paste(VCF.NAME ,"DBS.main.png", sep="."), 
#        width = 3000, height = 1200, units = "px", pointsize = 12)
    
# }, error = function(e)-999)

##################################################################################################
################################################################################################## 

# tryCatch({
#     
# p_dbs$p_dbs_cont
# ggsave(filename = paste(VCF.NAME ,"DBS.cont.png", sep="."), 
#        width = 3000, height = 1000, units = "px", pointsize = 12) 
# 
# }, error = function(e)-999)

# DBS
# Signature Contributions

################################################################################################### DBS
################################################################################################### 
################################################################################################### DBS
################################################################################################### Signature Contributions 

tryCatch({
    
sigs_dbs <-
  MutationalPatterns::get_known_signatures(muttype = "dbs") |>
  {
    \(sigs) sigrap::sig_contribution(mut_mat = dbs_counts, signatures = sigs)
  }()
    
sigs_dbs.df = as.data.frame(sigs_dbs)
sigs_dbs.df$sample = VCF.NAME
    
write.table(sigs_dbs.df, file=paste(VCF.NAME, "contributions.DBS.2020.txt", sep="."), 
sep="\t", col.names=TRUE,  row.names = FALSE, append = FALSE, quote = FALSE)
    
}, error = function(e)-999)
         
##################################################################################################
##################################################################################################
##################################################################################################
################################################################################################## 

tryCatch({

head(sigs_dbs.df, 2)

}, error = function(e)-999)

##################################################################################################
################################################################################################## HTML

tryCatch({

sigs_dbs |>
  sigrap::sig_contribution_table(type = "DBS", outdir = params$outdir) |>
  knitr::kable(format = "html") |>
  kableExtra::kable_styling(c("hover", "striped"), font_size = 12) |>
  kableExtra::scroll_box(height = "400px")
    
}, error = function(e)-999)
         
###################################################################################################
###################################################################################################
################################################################################################### INDELS
################################################################################################### PLOTS
###################################################################################################
################################################################################################### INDELS

tryCatch({

indel_counts <- sigrap::sig_count_indel(vcf_gr = gr, ref_genome = ref_genome)
p_indels <- sigrap::sig_plot_indel(indel_counts = indel_counts)

p_indels$p_indel_main / p_indels$p_indel_cont
ggsave(filename = paste(VCF.NAME ,"INDEL.png", sep="."), 
        width = 3200, height = 2200, units = "px", pointsize = 12)
    
     
}, error = function(e)-999)

##################################################################################################
##################################################################################################

# options(repr.plot.width = 14, repr.plot.height = 8)
# p_indels$p_indel_main 

# options(repr.plot.width = 14, repr.plot.height = 8)
# p_indels$p_indel_cont

# options(repr.plot.width = 14, repr.plot.height = 10)
# p_indels$p_indel_main / p_indels$p_indel_cont

##################################################################################################
################################################################################################## 

# tryCatch({

# p_indels$p_indel_main  
# ggsave(filename = paste(VCF.NAME ,"INDEL.main.png", sep="."), 
#       width = 3000, height = 2000, units = "px", pointsize = 12)
    
# p_indels$p_indel_cont
# ggsave(filename = paste(VCF.NAME ,"INDEL.cont.png", sep="."), 
#        width = 3000, height = 2000, units = "px", pointsize = 12)
#    
# }, error = function(e)-999)

################################################################################################### INDELS
################################################################################################### 
################################################################################################### SIGNATURE CONTRIBUTIONS
################################################################################################### INDELS 

tryCatch({

sigs_indel <-
  MutationalPatterns::get_known_signatures(muttype = "indel") |>
  {
    \(sigs) sigrap::sig_contribution(mut_mat = indel_counts, signatures = sigs)
  }()

    
sigs_indel.df = as.data.frame(sigs_indel)
sigs_indel.df$sample = VCF.NAME
 

write.table(sigs_indel.df, file=paste(VCF.NAME, "contributions.INDEL.2020.txt", sep="."), 
sep="\t", col.names=TRUE,  row.names = FALSE, append = FALSE, quote = FALSE)

    
}, error = function(e)-999)

##################################################################################################
##################################################################################################
##################################################################################################
################################################################################################## 

tryCatch({
head(sigs_indel.df, 2)   
    
}, error = function(e)-999)

################################################################################################### INDELS
################################################################################################### HTML 

tryCatch({
    
sigs_indel |>
  sigrap::sig_contribution_table(type = "ID", outdir = params$outdir) |>
  knitr::kable(format = "html") |>
  kableExtra::kable_styling(c("hover", "striped"), font_size = 12) |>
  kableExtra::scroll_box(height = "400px")
    
}, error = function(e)-999)

##################################################################################################
##################################################################################################

# to revisit the manual
# https://bioconductor.org/packages/release/bioc/html/MutationalPatterns.html

##################################################################################################
################################################################################################## 

# tryCatch({
#    
#    
# }, error = function(e)-999)

##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################