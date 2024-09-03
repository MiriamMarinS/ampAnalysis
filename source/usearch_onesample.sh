#!/bin/bash

# Amplicon analysis
# Miriam Marin Sanz, 2023

# Usearch pipeline


#Parameters: optimized from Step 1.
dif=58.94
pct=24.66
maxee=1.13
amp=$3 # minampsize
name_dir_usearch="/home/tonin/mimi/src"
project=$4

arg=$2 # name of the sample prefix
name_dir="./${project}/results_usearch" #Results directory name.

sample=$1 # path + sample prefix + sample suffix

# if log file exist, remove it
if [[ -f ${name_dir}/${arg}_log.txt ]]
then
   rm ${name_dir}/${arg}_log.txt
fi

#Usearch pipeline: merge.
echo "MERGING"
$name_dir_usearch/usearch9 -fastq_mergepairs $sample -relabel @ -fastq_maxdiffs $dif -fastq_maxdiffpct $pct -fastqout "${name_dir}/reads-${arg}.fq" &>> ${name_dir}/${arg}_log.txt
#Usearch pipeline: filter.
echo "FILTERING"
$name_dir_usearch/usearch9 -fastq_filter "${name_dir}/reads-${arg}.fq" -fastq_maxee $maxee -fastaout "${name_dir}/filtered-${arg}.fa" &>> ${name_dir}/${arg}_log.txt
#Usearch pipeline: uniques.
echo "UNIQUES"
$name_dir_usearch/usearch9 -fastx_uniques "${name_dir}/filtered-${arg}.fa" -fastaout "${name_dir}/${arg}_uniques.fa" -sizeout &>> ${name_dir}/${arg}_log.txt
#Usearch pipeline: unoise.
echo "UNOISE"
$name_dir_usearch/usearch9 -unoise2 "${name_dir}/${arg}_uniques.fa" -fastaout "${name_dir}/${arg}_denoised.fa" -minampsize $amp -log "${name_dir}/log-${arg}-unoise.txt" &>> ${name_dir}/${arg}_log.txt
#Usearch pipeline: searching.
echo "SEARCHING"
$name_dir_usearch/usearch9 -search_exact "${name_dir}/reads-${arg}.fq" -db "${name_dir}/${arg}_denoised.fa" -strand both -otutabout "${name_dir}/otutable_${arg}.txt" &>> ${name_dir}/${arg}_log.txt
