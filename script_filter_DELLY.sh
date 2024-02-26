#!/bin/bash

DELLY="/asclab-projects/spcg/working/wgsa-BT-resources/software/delly-0.7.5/delly_v0.7.5_CentOS5.4_x86_64bit"
SVPROPS="/asclab-projects/spcg/working/wgsa-BT-resources/software/svprops/src/svprops"
SAMPLEPROPS="/asclab-projects/spcg/working/wgsa-BT-resources/software/svprops/src/sampleprops"

REFERENCE_HG38_MALE="/asclab-projects/spcg/working/wgsa-BT-resources/genomes/genomeM.fa"
REFERENCE_HG38_FEMALE="/asclab-projects/spcg/working/wgsa-BT-resources/genomes/genomeF.fa"

LIST_CHRS_MALE="/asclab-projects/spcg/working/wgsa-BT-resources/genomes/hg38_male.interval_list"
LIST_CHRS_FEMALE="/asclab-projects/spcg/working/wgsa-BT-resources/genomes/hg38_female.interval_list"

## reading the INPUT files :

GENDER=""
FILE_with_SAMPLES_info=""

## an example :
## GENDER="MALE"
## FILE_with_SAMPLES_info="vcf.samples.for.delly.SPCG-AML449_5T.GRCh38p5M.MD.bam_vs_SPCG-AML449_8G.GRCh38p5M.MD.bam.tsv"
## SPCG-AML449_5T tumor
## SPCG-AML449_8G control

## an example :
#FILE_DEL="vcf.DELLY.SPCG-AML449_5T.GRCh38p5M.MD.bam_vs_SPCG-AML449_8G.GRCh38p5M.MD.bam.DEL.bcf"
#FILE_DUP="vcf.DELLY.SPCG-AML449_5T.GRCh38p5M.MD.bam_vs_SPCG-AML449_8G.GRCh38p5M.MD.bam.DUP.bcf"
#FILE_INS="vcf.DELLY.SPCG-AML449_5T.GRCh38p5M.MD.bam_vs_SPCG-AML449_8G.GRCh38p5M.MD.bam.INS.bcf"
#FILE_INV="vcf.DELLY.SPCG-AML449_5T.GRCh38p5M.MD.bam_vs_SPCG-AML449_8G.GRCh38p5M.MD.bam.INV.bcf"
#FILE_TRA="vcf.DELLY.SPCG-AML449_5T.GRCh38p5M.MD.bam_vs_SPCG-AML449_8G.GRCh38p5M.MD.bam.TRA.bcf"

FILE_DEL=""
FILE_DUP=""
FILE_INS=""
FILE_INV=""
FILE_TRA=""

###############################################
########### set up the reference genome

if [ "$GENDER" == "MALE" ]; then
   REFERENCE=$REFERENCE_HG38_MALE
fi

if [ "$GENDER" == "FEMALE" ]; then
   REFERENCE=$REFERENCE_HG38_FEMALE
fi

###################################################
###################################################
########### set up the list of chromosomes
###################################################

#if [ "$GENDER" == "MALE" ]; then
#   CHR=$LIST_CHRS_MALE
#fi

#if [ "$GENDER" == "FEMALE" ]; then
#   CHR=$LIST_CHRS_FEMALE
#fi

###################################################
###################################################
###################################################
########### filtering the VCF files : 

$DELLY filter \
-f somatic \
-t DEL \
-p \
-g $REFERENCE \
-m 15 \
-a 0.2 \
-r 0.75 \
-s $FILE_with_SAMPLES_info \
-o "${FILE_DEL%.bcf}.filtered.AF.0.2.bcf" \
$FILE_DEL

$DELLY filter \
-f somatic \
-t DUP \
-p \
-g $REFERENCE \
-m 400 \
-a 0.2 \
-r 0.75 \
-s $FILE_with_SAMPLES_info \
-o "${FILE_DUP%.bcf}.filtered.AF.0.2.bcf" \
$FILE_DUP

$DELLY filter \
-f somatic \
-t INS \
-p \
-g $REFERENCE \
-m 15 \
-a 0.2 \
-r 0.75 \
-s $FILE_with_SAMPLES_info \
-o "${FILE_INS%.bcf}.filtered.AF.0.2.bcf" \
$FILE_INS

$DELLY filter \
-f somatic \
-t INV \
-p \
-g $REFERENCE \
-m 400 \
-a 0.2 \
-r 0.75 \
-s $FILE_with_SAMPLES_info \
-o "${FILE_INV%.bcf}.filtered.AF.0.2.bcf" \
$FILE_INV

$DELLY filter \
-f somatic \
-t TRA \
-p \
-g $REFERENCE \
-m 0 \
-a 0.2 \
-r 0.75 \
-s $FILE_with_SAMPLES_info \
-o "${FILE_TRA%.bcf}.filtered.AF.0.2.bcf" \
$FILE_TRA

########################################################################
########################################################################
########################################################################

bcftools concat -a -O v \
-o "${FILE_with_SAMPLES_info%.tsv}.ALL.combined.AF.0.2.filtered.DEL.DUP.INS.INV.TRA.combined.vcf" \
"${FILE_DEL%.bcf}.filter.AF.0.2.bcf" \
"${FILE_DUP%.bcf}.filter.AF.0.2.bcf" \
"${FILE_INS%.bcf}.filter.AF.0.2.bcf" \
"${FILE_INV%.bcf}.filter.AF.0.2.bcf" \
"${FILE_TRA%.bcf}.filter.AF.0.2.bcf"


$SVPROPS "${FILE_with_SAMPLES_info%.tsv}.ALL.combined.AF.0.2.filtered.DEL.DUP.INS.INV.TRA.combined.vcf" > \
"${FILE_with_SAMPLES_info%.tsv}.ALL.combined.AF.0.2.filtered.DEL.DUP.INS.INV.TRA.combined.bedpe"

########################################################################
########################################################################
########################################################################

