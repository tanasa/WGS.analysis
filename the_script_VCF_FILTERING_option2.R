
library("VariantAnnotation")
library("VariantTools")
library("maftools")
library("GenVisR")
library("signeR")
library("SomaticSignatures")
library("TVTB")

###########################################################

vcf <- readVcf("AML.vcf",genome="hg38")

name_vcf <- basename("AML.vcf")

name_summary_file <- paste(trimws(name_vcf), "_results_summary", sep="")

###########################################################
###########################################################
###########################################################

# compressVcf <- bgzip("AML.vcf", "AML.vcf.bgzip")
# idx <- indexTabix(compressVcf, "vcf")
# tab <- TabixFile(compressVcf, idx)

####################################
####################################

# class(vcf)
# dim(vcf)

# header(vcf)
# geno(header(vcf))
# info(header(vcf))
# samples(header(vcf))

# info(vcf)
# info(header(vcf))

# geno(vcf)
# geno(header(vcf))

# rowRanges(vcf)
# fixed(vcf)

# fixed(vcf)
# ref(vcf)
# alt(vcf)
# qual(vcf)
# filt(vcf)

#########################

# colData(vcf)
# rowRanges(vcf)
# genome(vcf)
# seqlevels(vcf)

#######################

# isSNV(vcf)
# isInsertion(vcf)
# isDeletion(vcf)
# isTransition(vcf)

#######################

## res <- genotypeToSnpMatrix(vcf) ## Convert an array of genotype calls from 
## the "GT", "GP", "GL" or "PL" FORMAT field of a VCF file to a SnpMatrix.

######################################################################
######################################################################

#### Allele Fraction :

# head(geno(vcf)$AF[,"NORMAL"])
# head(geno(vcf)$AF[,"TUMOR"])
# head(geno(vcf)$AF[,c("NORMAL","TUMOR")])

#### Allele Depth :

# head(geno(vcf)$AD[,"NORMAL"])
# head(geno(vcf)$AD[,"TUMOR"])
# head(geno(vcf)$AD[,c("NORMAL","TUMOR")])

### to access the fields :
### head(geno(vcf)$AD[,"NORMAL"][1][[1]][1])
### head(geno(vcf)$AD[,"NORMAL"][1][[1]][2])
###################################
###################################
## $MIN_AD = 5;
## $MIN_AF = 0.05; 
## $MIN_RD = 10; 
## $MAX_RD = 100;
## $MAX_FS = 60;
## $MAX_SOR = 4;
###################################
###################################

######################################################################
######################################################################
######################################################################
################################### FILTERING a VCF OBJECT :


AF_AD_filters = function(x) {
                             af <- geno(x)$AF[,"TUMOR"] >= 0.05
                             ad <- min(as(geno(x)$AD[,"TUMOR"], "List")) >= 5
                             af & ad
                            }

AF_AD_rules <- FilterRules(list(AF_AD_filters = AF_AD_filters))

############################

FS_SOR_filters = function(x) {
                               fs  <- info(x)$FS <= 60
                               sor <- info(x)$SOR <= 4
                               fs & sor | !isSNV(x) 
                             }

FS_SOR_rules <- FilterRules(list(FS_SOR_filters = FS_SOR_filters))

##########################

vcf2 <- vcf[AF_AD_filters(vcf)]
vcf3 <- vcf2[FS_SOR_filters(vcf2)]

###############################

writeVcf(vcf3,"AML_AF_and_AD_FS_SOR_filtered_AS_OBJ_16june2017.vcf")


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

### printing the SUMMARIES :

###################################
################################### to count the total number of VARIANTS :

write.table("A summary on UN-FILTERED VCF file", file=name_summary_file, append=T, 
                          sep="\t", row.names=F, col.names=F)

count_variants <- dim(vcf)[1]
count_variants_file <- paste(c("total_number_variants :", count_variants), sep="\t")
write.table(count_variants_file, file=name_summary_file, append=T, sep="\t", 
                            row.names=F, col.names=F)


###################################
################################### to count the number of SNVS :

count_snp <- dim(vcf[isSNV(vcf)=="TRUE",])[1]
count_snp_file <- paste(c("total_number_snp :", count_snp), sep="\t") 
write.table(count_snp_file, file=name_summary_file, append=T, sep="\t", 
                            row.names=F, col.names=F)

###################################
################################### to count the number of INSERTIONS :

count_insertion <- dim(vcf[isInsertion(vcf)=="TRUE",])[1]
count_insertion_file <- paste(c("total_number_insertions :", count_insertion), sep="\t")
write.table(count_insertion_file, file=name_summary_file, append=T, sep="\t", 
                                  row.names=F, col.names=F)

###################################
################################### to count the number of DELETIONS :

count_deletion <- dim(vcf[isDeletion(vcf)=="TRUE",])[1]
count_deletion_file <- paste(c("total_number_deletions :", count_deletion), sep="\t")
write.table(count_deletion_file, file=name_summary_file, append=T, sep="\t", 
                                  row.names=F, col.names=F)


###########################################################################
###########################################################################
###########################################################################
###########################################################################

### to get the SUMMARY from a FILTERED VCF file :

vcf_filtered <- vcf3

###################################
################################### to count the total number of VARIANTS in the FILTERED FILE :

write.table("A summary on FILTERED VCF file", file=name_summary_file, append=T, sep="\t", 
                                              row.names=F, col.names=F)

count_variants_filtered <- dim(vcf_filtered)[1]
count_variants_filtered_file <- paste(c("total_number_variants :", count_variants_filtered), sep="\t")
write.table(count_variants_filtered_file, file=name_summary_file, append=T, sep="\t", 
                                          row.names=F, col.names=F)


###################################
################################### to count the number of SNVS :

count_snp_filtered <- dim(vcf_filtered[isSNV(vcf_filtered)=="TRUE",])[1]
count_snp_filtered_file <- paste(c("total_number_snp :", count_snp_filtered), sep="\t") 
write.table(count_snp_filtered_file, file=name_summary_file, append=T, sep="\t", 
                                     row.names=F, col.names=F)

###################################
################################### to count the number of INSERTIONS :

count_insertion_filtered <- dim(vcf_filtered[isInsertion(vcf_filtered)=="TRUE",])[1]
count_insertion_filtered_file <- paste(c("total_number_insertions :", count_insertion_filtered), sep="\t")
write.table(count_insertion_filtered_file, file=name_summary_file, append=T, sep="\t", 
                                           row.names=F, col.names=F)

###################################
################################### to count the number of DELETIONS :

count_deletion_filtered <- dim(vcf_filtered[isDeletion(vcf_filtered)=="TRUE",])[1]
count_deletion_filtered_file <- paste(c("total_number_deletions :", count_deletion_filtered), sep="\t")
write.table(count_deletion_filtered_file, file=name_summary_file, append=T, sep="\t", 
                                          row.names=F, col.names=F)


#######################################################################
#######################################################################

## Another way to filter a VCF OBJECT :

# vcf2_again <- subsetByFilter(vcf, filter=FilterRules(list(AF_AD_filter=AF_AD_filters)))
# vcf3_again <- subsetByFilter(vcf2_again, filter=FilterRules(list(FS_SOR_filtersr=FS_SOR_filters)))

# writeVcf(vcf3_again,"AML_out_AF_and_AD_FS_SOR_filtered_AS_OBJ_subsetbyfilter.vcf")

##########################################
##########################################

# readVcfAsVRanges(x, genome, param = ScanVcfParam(), ...): 
#                    Reads a VCF x directly into a VRanges; see readVcf for details on the arguments. 
#                    readVcfAsVRanges is an alternative syntax to as(readVcf(), "VRanges")

# NOTE: By default all INFO and FORMAT fields are read in with ScanVcfParam(). 
# The minimal information needed to create the VRanges can be specified as it follows:

# ScanVcfParam(fixed = "ALT", info = NA, geno = "AD")) 

# svp <- ScanVcfParam(geno=c("AD","AF"), samples="TUMOR")

# vcf2 <- readVcfAsVRanges("AML_expanded.vcf", "hg38", svp)

#######################################################################
#######################################################################