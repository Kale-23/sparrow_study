#!/bin/bash

#SBATCH --job-name=spades
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=256G
#SBATCH --array=1-390%20
#SBATCH --output=./output/spades/%j.out
#SBATCH --exclude=node\[117-118\]

# slurm batch script input
batch_list=/mnt/lz01/plachetzki/kcd1021/working_dir/913_project/spades_in_slurm

#run dir
file_in=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $batch_list)
file_in=$(realpath $file_in)
dir_in=$(dirname $(dirname $file_in))
dir_out=${dir_in}/mito_reads
mkdir -p $dir_out

echo $dir_out
echo $file_in
# actual command to run
./bam_to_reads.sh $file_in $dir_out
