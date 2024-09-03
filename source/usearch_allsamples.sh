#!/bin/bash

# Amplicon analysis
# Miriam Marin Sanz, 2023

# Usearch pipeline


#Parameters: optimized from Step 1.
dif=58.94
pct=24.66
maxee=1.13
amp=$2 # minampsize
name_dir_usearch="/home/tonin/mimi/src"
project=$3

name_dir="./${project}/results_usearch" #Results directory name.

input_files=$1 # path

# if log file exist, remove it
if [[ -f ${name_dir}/${project}_log.txt ]]
then
   rm ${name_dir}/${project}_log.txt
fi

#Usearch pipeline: merge.
echo "MERGING"
$name_dir_usearch/usearch9 -fastq_mergepairs $input_files/*R1*.fastq -relabel @ -fastq_maxdiffs $dif -fastq_maxdiffpct $pct -fastqout "${name_dir}/reads-${project}.fq" &>> ${name_dir}/${project}_log.txt
#Usearch pipeline: filter.
echo "FILTERING"
$name_dir_usearch/usearch9 -fastq_filter "${name_dir}/reads-${project}.fq" -fastq_maxee $maxee -fastaout "${name_dir}/filtered-${project}.fa" &>> ${name_dir}/${project}_log.txt
#Usearch pipeline: uniques.
echo "UNIQUES"
$name_dir_usearch/usearch9 -fastx_uniques "${name_dir}/filtered-${project}.fa" -fastaout "${name_dir}/${project}_uniques.fa" -sizeout &>> ${name_dir}/${project}_log.txt
#Usearch pipeline: unoise.
echo "UNOISE"
$name_dir_usearch/usearch9 -unoise2 "${name_dir}/${project}_uniques.fa" -fastaout "${name_dir}/${project}_denoised.fa" -minampsize $amp -log "${name_dir}/log-${project}-unoise.txt" &>> ${name_dir}/${project}_log.txt
#Usearch pipeline: searching.
echo "SEARCHING"
$name_dir_usearch/usearch9 -search_exact "${name_dir}/reads-${project}.fq" -db "${name_dir}/${project}_denoised.fa" -strand both -otutabout "${name_dir}/otutable_${project}.txt" &>> ${name_dir}/${project}_log.txt
