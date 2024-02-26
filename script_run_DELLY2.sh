#!/bin/bash

DELLY="/N/users/btanasa/software/delly-0.7.5/delly_v0.7.5_CentOS5.4_x86_64bit" 
SAMTOOLS="samtools"

DBSNP138="/N/projects/spcg/subjects/resources/Homo_sapiens_assembly38.dbsnp138.vcf"
COSMIC76="/N/projects/spcg/subjects/resources/COSMIC76-coding-and-noncoding.21march2016.vcf"

REFERENCE_HG38_MALE="/N/projects/spcg/subjects/resources/genomeM.fa"
REFERENCE_HG38_FEMALE="/N/projects/spcg/subjects/resources/genomeF.fa"

LIST_CHRS_MALE="/N/projects/spcg/subjects/resources/hg38_male.interval_list"
LIST_CHRS_FEMALE="/N/projects/spcg/subjects/resources/hg38_female.interval_list"

GENDER=$1
NORMAL_MD=$2
TUMOR_MD=$3

if [ "$GENDER" == "MALE" ]; then  
   REFERENCE=$REFERENCE_HG38_MALE
fi

if [ "$GENDER" == "FEMALE" ]; then 
   REFERENCE=$REFERENCE_HG38_FEMALE
fi 

$DELLY call -t DEL -o vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.DEL.bcf -g $REFERENCE $TUMOR_MD $NORMAL_MD 
$DELLY call -t INS -o vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.INS.bcf -g $REFERENCE $TUMOR_MD $NORMAL_MD 
$DELLY call -t INV -o vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.INV.bcf -g $REFERENCE $TUMOR_MD $NORMAL_MD 
$DELLY call -t DUP -o vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.DUP.bcf -g $REFERENCE $TUMOR_MD $NORMAL_MD 
$DELLY call -t TRA -o vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.TRA.bcf -g $REFERENCE $TUMOR_MD $NORMAL_MD 

bcftools view vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.DEL.bcf > vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.DEL.vcf
bcftools view vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.INS.bcf > vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.INS.vcf
bcftools view vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.INV.bcf > vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.INV.vcf
bcftools view vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.DUP.bcf > vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.DUP.vcf
bcftools view vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.TRA.bcf > vcf.DELLY.${TUMOR_MD}_vs_${NORMAL_MD}.TRA.vcf



