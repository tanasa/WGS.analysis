# CHORD
# https://github.com/UMCUGenetics/CHORD

library(BSgenome) 
library(BSgenome.Hsapiens.UCSC.hg38)
library(randomForest)
library(VariantAnnotation)
library(xtable)
library(kableExtra)

# randomForest is required by CHORD
# snvindel_vcf =
# sv_vcf =
# sample = 

# SVINDEL = "strelka2.consensus.snpEff.annotated.intersected.filtered.hg38.vcf" 
# SV = "manta.somatic.svs.hg38.vcf"
# SAMPLE ="EE"

args = commandArgs(trailingOnly=TRUE)

SVINDEL= args[1]
SV = args[2]
SAMPLE = args[3]

SAMPLE = paste(SAMPLE, "chord", sep=".")

# to test whether the files have been read correctly

svindel <- readVcf(SVINDEL, "hg38")
header(svindel)
head(rowRanges(svindel), 3)

sv <- readVcf(SV, "hg38")
header(sv)
head(rowRanges(sv), 3)

snvindel_vcf <- SVINDEL
sv_vcf <- SV

sv_df <- sigrap::chord_mantavcf2df(sv_vcf) # prepare SV VCF as data.frame
# dir.create(sample)

res <- sigrap::chord_run(
  vcf.snv = snvindel_vcf,
  df.sv = sv_df,
  sv.caller = "manta",
  # vcf.sv = sv_vcf,  ######## alternative
  ref.genome = "hg38",
  sample.name = SAMPLE,
  verbose = TRUE
)

x = data.frame(res$prediction)

# kbl(x)
# y = cbind(var = colnames(res$prediction), value = unlist(res$prediction[1, ])) |>
#    dplyr::as_tibble(.name_repair = "check_unique") |>
#    knitr::kable(caption = paste(sample, "CHORD results", sep="."))
    
y = cbind(var = colnames(res$prediction), value = unlist(res$prediction[1, ])) |>
     dplyr::as_tibble(.name_repair = "check_unique") 

# y %>% kable_paper("hover", full_width = F)
# y %>% kable_classic(full_width = F, html_font = "Cambria")

write.table(as.data.frame(y), file=paste(SAMPLE, "STRELKA.results.CHORD.txt", sep="."),
                              quote=FALSE, sep="\t", col.names=FALSE)
