#!/bin/bash

## STRELKA="/N/users/btanasa/Software/strelka_workflow-1.0.14/bin/configureStrelkaWorkflow.pl"
## CONFIG="/N/users/btanasa/Software/strelka_workflow-1.0.14/bin/strelka_config_bwa_default.ini"

MANTA="/N/users/btanasa/software/manta-0.29.6.centos5_x86_64/bin/configManta.py"
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


## GENDER="MALE"
## NORMAL_MD="SPCG-WT436_5G.GRCh38p5M.MD.bam"
## TUMOR_MD="SPCG-WT436_6T.GRCh38p5M.MD.bam"
## TUMOR_MD="SPCG-WT436_10T.GRCh38p5M.MD.bam"


$MANTA \
--normalBam=$NORMAL_MD \
--tumorBam=$TUMOR_MD \
--referenceFasta=$REFERENCE \
--runDir=./vcf-manta-0296-bwa_${GENDER}_${TUMOR_MD}_vs_${NORMAL_MD}

cd ./vcf-manta-0296-bwa_${GENDER}_${TUMOR_MD}_vs_${NORMAL_MD}
runWorkflow.py -m local -j 8


