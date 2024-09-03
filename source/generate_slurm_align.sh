#!/bin/bash

project_name=$1
### Directory containing raw FASTQ files to process
data_dir=$2

### Directory to send non-log outputs to
### (If specified dir does not exist, script creates it)
out_dir="./${project_name}/indels"

### Directory to send slurm logs to
### (If specified dir does not exist, script creates it)
log_dir="./${project_name}/indels/logs/"

### Partition(s) to submit jobs to
### (If want to specify >1 partition then comma separate them with no spaces)
partitions=""

### Max. memory allowed per array task
### (This is split across all CPUs allocated to that array task)
memory="24G"

### CPUs per array task
num_CPUs="1"

### Time per task (days-hours:minutes:seconds)
permitted_time="0-12:00"

### Name you want to give to SLURM master job
job_name="${project_name}"

PATHSCRIPTS="/home/tonin/mimi/protoplasts_alpha/Scripts/ampAnalysisv2/source"
##### MAIN BODY #####
#####################

### Variable set to current date and time so the SLURM array script generated
### always has a unique ID
script_uniq_ID=$(date +"%Y-%m-%d_%H_%M")

### Save a list of all files in specified directory into file_list.tmp.txt
### (Could pipe output from 'ls' into 'grep' to filter for specific file types)
ls -1 $data_dir | grep denoised > ${out_dir}/file_list_align.tmp.txt

array_length=$(($(wc -l < ${out_dir}/file_list_align.tmp.txt) - 1))

### Create directory for main outputs
### If it already exists will get an error, I use "2> /dev/null" to delete this
### message because it is not useful information; error or not we still get the
### desired output
mkdir $out_dir 2> /dev/null

### Create directory for slurm logs
### If it already exists will get an error, "2> /dev/null" deletes this message
mkdir $log_dir 2> /dev/null

### Create the SLURM array script and print SBATCH messages to it
### SBATCH parameters are influenced by the user specified inputs above

echo "#!/bin/bash
#SBATCH --job-name="$job_name"
#SBATCH --array=0-"$array_length"
#SBATCH --output=$log_dir/"$job_name"_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task="$num_CPUs"
#SBATCH --mem-per-cpu=12G
#SBATCH --time="$permitted_time"
" > ${out_dir}/array_script_align.slurm

### SCRIPT

echo "
mapfile -t FILE_ARRAY < ${out_dir}/file_list_align.tmp.txt

sample_mate1_filename=\${FILE_ARRAY[ \$SLURM_ARRAY_TASK_ID ]}
sample_prefix=\${sample_mate1_filename%_denoised.fa}
echo 
bwa mem -t 18 -M $3 ${data_dir}\${sample_prefix}_denoised.fa > ${out_dir}/\${sample_prefix}_bwa.sam
mkdir $out_dir/temp 2> /dev/null
mkdir $out_dir/temp/\${sample_prefix} 2> /dev/null # Create a folder for tmp files
bbmap.sh -Xmx4g in=${data_dir}\${sample_prefix}_denoised.fa out=${out_dir}/\${sample_prefix}_bbmap.sam sam=1.3 ref=$4 path=${out_dir}/temp/\${sample_prefix}      
echo ""
echo Finished this sub-task!
" >> ${out_dir}/array_script_align.slurm