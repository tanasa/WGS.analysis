
#!/usr/bin/bash` 

ANNOVAR="/home/bogdan/ANNOVAR/annovar/annotate_variation.pl"
ANNOVAR_DATABASES="/home/bogdan/ANNOVAR/annovar/databases/"

SAMPLE=$1

awk 'FNR>1 || NR==1' $SAMPLE | grep -v "seqnames" | awk '{print $1"\t"$2"\t"$3"\t"0"\t"0"\t"$9}' > "${SAMPLE}.for.annovar"

$ANNOVAR "${SAMPLE}.for.annovar" $ANNOVAR_DATABASES --buildver hg38 --separate --expandbin 10000000