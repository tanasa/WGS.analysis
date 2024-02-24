#!/bin/bash

# ls -1 *vcf | grep "manta" | cut -d '.' -f 1 > the.list.of.ID.txt

# to have all the strings in an array
# IFS=$'\n' read -d '' -r -a arr < the.list.of.ID.txt

# we can read the strings in an array in two ways, 
# either with mapfile  or with readarray

# mapfile -t arr < the.list.of.ID.txt
# echo "${arr[@]}"
# echo "${arr[1]}"
# readarray -t arr2 < the.list.of.ID.txt
# echo "${arr2[1]}"
# echo "${arr2[@]}"

# to list the strings that are collected in the array
# for index in "${!arr[@]}";
# do
# echo "$index -> ${arr[$index]}"
# done

# echo "${arr[@]}"

# another way to list the elements in the array
# for (( i=0; i<${#arr[@]}; i++ ));
# do
# echo ${arr[$i]}
# done

# another way to list the elements in the array

# for index in "${!arr[@]}";
# do

# echo "${arr[$index]}"
# sample="${arr[$index]}"
# snv="${sample}.tnscope.eff_kgsnp_kgindel_mgindel_cosmic_exac_dbsnp_gnomad.intersected.hg38.vcf"
# sv="${sample}.manta.tumor.svs.hg38.vcf"

# echo $sample
# echo $snv
# echo $sv

# Rscript --vanilla the_script_analysis_CHORD.R $snv $sv $sample

# done

# another  way to generate the data

file="the.list.of.ID.txt"
lines=`cat $file`

for line in $lines; do

echo "$line"
sample="$line"
snv="${line}.somatic.strelka2.only.mutations.hg38.vcf"
sv="${line}.manta.tumor.svs.hg38.vcf"

echo $sample
echo $snv
echo $sv

Rscript --vanilla the_script_analysis_CHORD_STRELKA.R $snv $sv $sample

done
