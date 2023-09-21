###########################################################################
###########################################################################
###########################################################################
###########################################################################

library("VariantAnnotation")

# library("GenVisR")
# library("VariantTools")
# library("maftools")
# library("SomaticSignatures")
# library("signeR")

###########################################################

vcf <- readVcf("AML.vcf",genome="hg38")
name_vcf <- basename("AML.vcf")

# to trim the WHITE SPACES :
# returns string w/o leading whitespace
# trim.leading <- function (x)  sub("^\\s+", "", x)

# returns string w/o trailing whitespace
# trim.trailing <- function (x) sub("\\s+$", "", x)

# returns string w/o leading or trailing whitespace
# trim <- function (x) gsub("^\\s+|\\s+$", "", x)

# trim <- function (x) gsub("^\\s+|\\s+$", "", x)
# name_vcf <- trim(name_vcf)

name_summary_file <- paste(trimws(name_vcf), "_results_summary", sep="")

##############################################

compressVcf <- bgzip("AML.vcf", "AML.vcf.bgzip")
idx <- indexTabix(compressVcf, "vcf")
tab <- TabixFile(compressVcf, idx)

####################################
####################################

# class(vcf)
# dim(vcf)  ## it prints the number of variants

# header(vcf)
# geno(header(vcf))
# info(header(vcf))
# samples(header(vcf))

# rowRanges(vcf)                     ## ROW_RANGES FIELD 

# info(vcf)                          ## INFO FIELD
# head(info(vcf))

# geno(vcf)                          ## GENO FIELD  

# fixed(vcf)                         ## FIXED FIELD  

# ref(vcf)                           ## REF FIELD 
# alt(vcf)                           ## ALT FIELD 
# qual(vcf)                          ## QUAL FIELD
# filt(vcf)                          ## FILTER FIELD

#########################

# colData(vcf)                       ## samples in the GENO FIELD 
# rowRanges(vcf)                     ## listing the ROWS 
# genome(vcf)
# seqlevels(vcf)

#######################

# isSNV(vcf)                         $
# isInsertion(vcf)
# isDeletion(vcf)

# isTransition(vcf)
# isSubstitution(vcf)

#######################

## res <- genotypeToSnpMatrix(vcf) # Convert an array of genotype calls  
## the "GT", "GP", "GL" or "PL" FORMAT field of a VCF file to a SnpMatrix.

######################################################################
######################################################################

#### Allele Fraction :

# head(geno(vcf)$AF[,"NORMAL"])
# head(geno(vcf)$AF[,"TUMOR"])
# head(geno(vcf)$AF[,c("NORMAL","TUMOR")])

##### Allele Depth :

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
################################### to count the total number of VARIANTS :

write.table("A summary on UN-FILTERED VCF file", file=name_summary_file, append=T, sep="\t", 
                            row.names=F, col.names=F)

count_variants <- dim(vcf)[1]
count_variants_file <- paste(c("total_number_variants :", count_variants), sep="\t")
write.table(count_variants_file, file=name_summary_file, append=T, sep="\t", 
                            row.names=F, col.names=F)

###################################
################################### to count the number of SNPS :

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

###################################
#########################################################################
#########################################################################
################################### FILTERING the VCF FILE :

AF_AD_filters = function(x) {
                             af <- geno(x)$AF[,"TUMOR"] >= 0.05
                             ad <- min(as(geno(x)$AD[,"TUMOR"], "List")) >= 5
                             af & ad
                            }

AF_AD_rules <- FilterRules(list(AF_AD_filters = AF_AD_filters))

#########################################################################
#########################################################################
##### To filter the VCF file, use the filter as the argument to filterVcf() :

vcf_filtered <- filterVcf( "AML.vcf.bgzip", 
                           "hg38", 
                           "AML_out_AF_and_AD_filtered.vcf", 
                            filters=AF_AD_rules)

#######################################################################
#######################################################################

compressVcf <- bgzip("AML_out_AF_and_AD_filtered.vcf", 
                     "AML_out_AF_and_AD_filtered.vcf.bgzip")

idx <- indexTabix(compressVcf, "vcf")
tab <- TabixFile(compressVcf, idx)

#######################################################################
#######################################################################

FS_SOR_filters = function(x) {
                               fs  <- info(x)$FS <= 60
                               sor <- info(x)$SOR <= 4
                               fs & sor | !isSNV(x) 
                             }

FS_SOR_rules <- FilterRules(list(FS_SOR_filters = FS_SOR_filters))

vcf_filtered <- filterVcf( "AML_out_AF_and_AD_filtered.vcf.bgzip", 
                           "hg38", 
                           "AML_out_AF_and_AD_FS_and_SOR_filtered_way1.vcf", 
                            filters=FS_SOR_rules)

###########################################################################
###########################################################################

### to get the SUMMARY from a FILTERED VCF file, we shall do :

vcf_filtered <- readVcf("AML_out_AF_and_AD_FS_and_SOR_filtered_way1.vcf",genome="hg38")

###########################################################################
###########################################################################
###################################
################################### to count the total number of VARIANTS in the FILTERED FILE :

write.table("A summary on FILTERED VCF file", file=name_summary_file, append=T, sep="\t", 
                                              row.names=F, col.names=F)

count_variants_filtered <- dim(vcf_filtered)[1]
count_variants_filtered_file <- paste(c("total_number_variants :", count_variants_filtered), sep="\t")
write.table(count_variants_filtered_file, file=name_summary_file, append=T, sep="\t", 
                                          row.names=F, col.names=F)

###########################################################################
###########################################################################
###################################
################################### to count the number of SNVS:

count_snp_filtered <- dim(vcf_filtered[isSNV(vcf_filtered)=="TRUE",])[1]
count_snp_filtered_file <- paste(c("total_number_snp :", count_snp_filtered), sep="\t") 
write.table(count_snp_filtered_file, file=name_summary_file, append=T, sep="\t", 
                                     row.names=F, col.names=F)

###########################################################################
###########################################################################
###################################
################################### to count the number of INSERTIONS :

count_insertion_filtered <- dim(vcf_filtered[isInsertion(vcf_filtered)=="TRUE",])[1]
count_insertion_filtered_file <- paste(c("total_number_insertions :", count_insertion_filtered), sep="\t")
write.table(count_insertion_filtered_file, file=name_summary_file, append=T, sep="\t", 
                                           row.names=F, col.names=F)

###########################################################################
###########################################################################
###################################
################################### to count the number of DELETIONS :

count_deletion_filtered <- dim(vcf_filtered[isDeletion(vcf_filtered)=="TRUE",])[1]
count_deletion_filtered_file <- paste(c("total_number_deletions :", count_deletion_filtered), sep="\t")
write.table(count_deletion_filtered_file, file=name_summary_file, append=T, sep="\t", 
                                          row.names=F, col.names=F)


###########################################################################
###########################################################################
###########################################################################
###########################################################################
## To filter a VCF object, use subsetByFilter(vcf, filter) and see next ;)


# readVcfAsVRanges(x, genome, param = ScanVcfParam(), ...): 
#                     Reads a VCF x directly into a VRanges; see readVcf for details on the arguments. 
#                     readVcfAsVRanges is an alternative syntax to as(readVcf(), "VRanges")

# NOTE: By default all INFO and FORMAT fields are read in with ScanVcfParam(). 
# The minimal information needed to create the VRanges can be specified as it follows:

# ScanVcfParam(fixed = "ALT", info = NA, geno = "AD")) 

# svp <- ScanVcfParam(geno=c("AD","AF"), samples="TUMOR")


# vcf2 <- readVcfAsVRanges("AML_expanded.vcf", "hg38", svp)

###########################################################################
###########################################################################
###########################################################################
###########################################################################