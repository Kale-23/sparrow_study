#!/bin/bash
#SBATCH --job-name=stuff
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --array=1-390%15
#SBATCH --output=./output/mito/%j.out
#SBATCH --exclude=node\[117-118\]

#echo "Usage: $0 <RUN_DIR> <THREADS> <FILE> <REF_INDEX> [READS]"
#echo "  RUN_DIR       : Directory of output" echo "  THREADS       : Number of threads to be used"
#echo "  FILE          : Trimmomatic file to process (excluding _R[1-2].fastq.gz)"
#echo "  REF_INDEX     : Location of Reference file"
#echo "  READS      : Number of reads to direct seqtk to keep. Default is all"
#echo "  -h        Display this help message"

# slurm batch script input
batch_list=/mnt/lz01/plachetzki/kcd1021/working_dir/913_project/Mmelodia_slurm_input_with_reads

#run dir
cov_dir=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $batch_list)
run_dir=cov_"$cov_dir"

#threads
threads=24

#file
file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $batch_list)

#ref index
#ref_index=/mnt/lz01/plachetzki/kcd1021/working_dir/913_project/ref_indexes/GCF_035770615.1_bMelMel2.pri_genomic.fna
ref_index=/mnt/lz01/plachetzki/kcd1021/working_dir/913_project/ref_indexes/just_mito/CM069588.1.fna

#reads
reads=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $batch_list)

#run script
./full_project_script.sh $run_dir $threads $file $ref_index $reads

#### just for running flagstat ###
#FILE=$sample
#HOME_DIR="/mnt/lz01/plachetzki/kcd1021/working_dir/913_project/$ouput_dir"
#ALIGN_DIR=$HOME_DIR/alignment_out
#OUTPUT_SORTED_BAM=$ALIGN_DIR/"$FILE"aligned_reads_sorted.bam
#
#FLAGSTAT_DIR=$HOME_DIR/flagstat
#mkdir -p $FLAGSTAT_DIR
#OUTPUT_FLAGSTAT=$FLAGSTAT_DIR/"$FILE"flagstat.tsv
#
#samtools flagstat --output-fmt tsv --threads 24 $OUTPUT_SORTED_BAM >$OUTPUT_FLAGSTAT
