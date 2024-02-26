
#!/bin/bash

## the bash script reads the gender (MALE, FEMALE), the NORMAL and the TUMOR files (BAM)
## it runs after LUMPY

PICARD="java -jar /N/users/btanasa/Software/picard-tools-1.140/picard.jar"
SAMTOOLS="samtools"
GATK="java -Xmx8g -jar /N/users/btanasa/GATK-3.5/GenomeAnalysisTK.jar"
MUTECT1="java -Xmx8g -jar /N/users/btanasa/Software/mutect-1.1.7/mutect-1.1.7.jar" 

DBSNP138="/N/projects/spcg/subjects/resources/Homo_sapiens_assembly38.dbsnp138.vcf"
COSMIC76="/N/projects/spcg/subjects/resources/COSMIC76-coding-and-noncoding.21march2016.vcf"

REFERENCE_HG38_MALE="/N/projects/spcg/subjects/resources/genomeM.fa"
REFERENCE_HG38_FEMALE="/N/projects/spcg/subjects/resources/genomeF.fa"

SVTYPER="/N/users/btanasa/software/svtyper-master_23nov2016/svtyper"

VCF=$1
NORMAL=$2
TUMOR=$3 

if [ "$GENDER" == "MALE" ]; then  
   REFERENCE=$REFERENCE_HG38_MALE
fi

if [ "$GENDER" == "FEMALE" ]; then 
   REFERENCE=$REFERENCE_HG38_FEMALE
fi 

echo "we do run the analysis with SVTYPER :"


$SVTYPER -i $VCF -B $TUMOR,$NORMAL -l "${VCF}.json" > "${VCF}.with.SVTYPER.vcf" 

