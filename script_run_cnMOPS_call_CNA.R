#!/usr/bin/env Rscript

library(cn.mops)
library(parallel)
library(snow)
library(RUnit)
library(BiocParallel)

#############################################################
#############################################################

## a change : 14feb2018
## removing <mode="paired"> from the function getReadCountsFromBam()

# registered()
# register(MulticoreParam(workers = 1))

####################################################################################################################
####################################################################################################################

args <- commandArgs(TRUE)

## args=(commandArgs(TRUE))
## args is now a list of character vectors

## TUMOR <- args[1]         
## GERMLINE <- args[2]
## CHR <- args[3] 

## the order of the files is TUMOR, and GERMLINE :

TUMOR <- ""     
GERMLINE <- ""
CHR <- args[1] 

## First check to see if arguments are passed.
## if(length(args)==0) { 
##                     stop("no args specified") 
##                    }
## Then cycle through each element of the list and evaluate the expressions.
## for(i in 1:length(args)){
##                         print(args[[i]])
##                         eval(parse(text=args[[i]]))
##                        }
## print(TUMOR)
## print(GERMLINE)
## print(CHR)

####################################################################################################################
####################################################################################################################

tumor <- getReadCountsFromBAM(TUMOR, refSeqName=CHR, WL=10000, parallel=12)
normal <- getReadCountsFromBAM(GERMLINE, refSeqName=CHR, WL=10000, parallel=12)

### a special normalization i.e. POISSON because the tumor has made large CNVs

X <- tumor
values(X) <- cbind(values(tumor),values(normal)) 

### NORMALIZATION :

X.mode <- normalizeGenome(X, normType="poisson")

### Parameter settings for tumor: 
### - norm=0, because we already have normalized
### - integer copy numbers higher than 8 allowed
### - DNAcopy as segmentation algorithm.

### the position 1 is for TUMOR, the position 2 is for NORMAL.

ref_analysis_norm0 <- referencecn.mops(X.mode[,1], X.mode[,2],
	norm=0, 
	I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 8, 16, 32, 64), 
	classes = c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN7","CN8","CN16","CN32","CN64","CN128"),
	segAlgorithm="DNAcopy")

####################################################################################################################
######################################################################## if there are results of the CNV analysis :
####################################################################################################################

name_file_results_analysis <- paste("z.cn.mops.segments.normalization.POISSON", basename(TUMOR), basename(GERMLINE), basename(CHR), "analysis.results", sep=".")

if (length(cnvs(ref_analysis_norm0))==0) {
     writeLines("no CNV detected",name_file_results_analysis)
} else {
     writeLines("CNV were detected",name_file_results_analysis)
}

####################################################################################################################
##########################################################################
####################################################################################################################

resCNMOPS <- calcIntegerCopyNumbers(ref_analysis_norm0)
resCNMOPS <- cn.mops:::.replaceNames(resCNMOPS, basename(TUMOR) ,"TUMOR")
 
###########################################################################
## name_file <- paste(args[1], args[2], args[3], ".png", sep="_")
## name_file_display <- paste("z.cn.mops.seqplot.", basename(TUMOR), ".", basename(GERMLINE), ".", basename(CHR), ".png")
###########################################################################

name_file_display <- paste("cn.mops.seqplot", basename(TUMOR), basename(GERMLINE), basename(CHR), "png", sep=".")

png(name_file_display)
### segplot(resCNMOPS)
segplot(resCNMOPS, seqnames=CHR, zlcol="blue", segcol="red")
dev.off()

####################################################################################################################
###########################################################################
####################################################################################################################

### cnvs(resCNMOPS) # look at individual CNV regions
### cnvr(resCNMOPS) # look at REGIONS 

####################################################################################################################
############### SAVING THE RESULTS :
####################################################################################################################

results.segm  <- as.data.frame(segmentation(resCNMOPS))
results.CNVs  <- as.data.frame(cnvs(resCNMOPS))
results.CNVRegions  <- as.data.frame(cnvr(resCNMOPS))

####################################################################################################################
############### NAMING THE FILES :
####################################################################################################################

name_file_results_segm <- paste("cn.mops.segments.normalization.POISSON", basename(TUMOR), basename(GERMLINE), basename(CHR), "txt", sep=".")
name_file_results.CNVs <- paste("cn.mops.cnvs_shorts.normalization.POISSON", basename(TUMOR), basename(GERMLINE), basename(CHR), "txt", sep=".")
name_file_results.CNVRegions <- paste("cn.mops.cnvs_regions.normalization.POISSON", basename(TUMOR), basename(GERMLINE), basename(CHR), "txt", sep=".") 

####################################################################################################################
############### WRITING THE TABLES :
####################################################################################################################

write.table(results.segm, file=name_file_results_segm, sep="\t", row.names=F)
write.table(results.CNVs, file=name_file_results.CNVs, sep="\t", row.names=F)
write.table(results.CNVRegions, file=name_file_results.CNVRegions, sep="\t", row.names=F) 

##############################################################################################
##############################################################################################

### a series of functions :
### gr(ref_analysis_norm0)

### seqnames(gr(ref_analysis_norm0))
### ranges(gr(ref_analysis_norm0)) 
### strand(gr(ref_analysis_norm0))
### elementMetadata(gr(ref_analysis_norm0))
### seqinfo(gr(ref_analysis_norm0)) 
### metadata(gr(ref_analysis_norm0)) 

### normalizedData(ref_analysis_norm0)
### localAssessments(ref_analysis_norm0)
### individualCall(ref_analysis_norm0) 
### posteriorProbs(ref_analysis_norm0)

### cnvs(ref_analysis_norm0) ## : length(cnvs(ref_analysis_norm0)) == 0
### cnvr(ref_analysis_norm0) ## : length(cnvr(ref_analysis_norm0)) == 0

### segmentation(ref_analysis_norm0)

### integerCopyNumber(ref_analysis_norm0) # for a chromosome without CNVm these appear as CN2

### params(ref_analysis_norm0) # 

##############################################################################################
##############################################################################################