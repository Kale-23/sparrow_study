#!/bin/bash

usage() {
	echo "Usage: $0 <RUN_DIR> <THREADS> <FILE> <REF_INDEX> [READS]"
	echo "  RUN_DIR       : Directory of output"
	echo "  THREADS       : Minimum number of threads to be used"
	echo "  FILE          : Trimmomatic file to process (excluding _R[1-2].fastq.gz)"
	echo "  REF_INDEX     : Location of Reference file"
	echo "  READS      : Number of reads to direct seqtk to keep. Default is all"
	echo "  -h        Display this help message"
	exit 1
}

echo starting at $(date)
# useage message if no arguments or if -h
if [ $# -lt 3 ] || [ "$1" == "-h" ]; then
	usage
fi

RUN_DIR="$1"
THREADS="$2"
FILE="$3"
REF_INDEX="$4"

if [[ $5 == "" ]]; then
	READS="full"
else
	READS=$5
fi

echo processing $FILE

# used in all parts
#REF_DIR="/mnt/lz01/plachetzki/kcd1021/working_dir/913_project/ref_indexes"
#HOME_DIR="/mnt/lz01/plachetzki/kcd1021/working_dir/913_project/subs/$RUN_DIR"
HOME_DIR="/mnt/lz01/plachetzki/kcd1021/working_dir/913_project/mito/$RUN_DIR"
mkdir -p $HOME_DIR
TRIM_DIR="/mnt/lz01/plachetzki/kcd1021/working_dir/913_project/trimmomatic_files"

module load linuxbrew/colsa
if [[ $READS != "full" ]]; then
	echo subsampling reads using seqtk
	SUBSAMPLE_READS_OUT=$HOME_DIR/subsampled_reads
	mkdir -p $SUBSAMPLE_READS_OUT

	SEED=$(($RANDOM % 256))
	echo subsampling at $(basename $HOME_DIR) with $READS
	echo seed: $SEED
	SUB_READ_1=$SUBSAMPLE_READS_OUT/"$FILE"R1.fastq
	SUB_READ_2=$SUBSAMPLE_READS_OUT/"$FILE"R2.fastq
	seqtk sample -s$SEED $TRIM_DIR/"$FILE"R1.fastq.gz $READS >$SUB_READ_1
	seqtk sample -s$SEED $TRIM_DIR/"$FILE"R2.fastq.gz $READS >$SUB_READ_2

	echo zipping trimmed fastqs $(date)
	gzip $SUB_READ_1
	gzip $SUB_READ_2

	ALIGN_READS_IN=$SUBSAMPLE_READS_OUT
else
	ALIGN_READS_IN=$TRIM_DIR
fi

ALIGN_DIR=$HOME_DIR/alignment_out
mkdir -p $ALIGN_DIR
ALIGN_LOG_DIR=$HOME_DIR/log_files/alignment_logs
mkdir -p $ALIGN_LOG_DIR

READ1="$ALIGN_READS_IN"/"$FILE"R1.fastq.gz
READ2="$ALIGN_READS_IN"/"$FILE"R2.fastq.gz

OUTPUT_LOG=$ALIGN_LOG_DIR/"$FILE"aligned.log
OUTPUT_SAM=$ALIGN_DIR/"$FILE"aligned_reads.sam
OUTPUT_BAM=$ALIGN_DIR/"$FILE"aligned_reads.bam
OUTPUT_SORTED_BAM=$ALIGN_DIR/"$FILE"aligned_reads_sorted.bam

echo running bwa mem $(date)
bwa mem -t 24 $REF_INDEX $READ1 $READ2 >$OUTPUT_SAM 2>$OUTPUT_LOG

echo running samtools view $(date)
samtools view --threads 24 -b -o $OUTPUT_BAM $OUTPUT_SAM &&
	rm $OUTPUT_SAM

echo running samtools sort $(date)
samtools sort --threads 24 -l 9 -o $OUTPUT_SORTED_BAM $OUTPUT_BAM &&
	rm $OUTPUT_BAM

echo running samtools index $(date)
samtools index -@ 24 $OUTPUT_SORTED_BAM

FLAGSTAT_DIR=$HOME_DIR/flagstat
mkdir -p $FLAGSTAT_DIR

echo running samtools flagstat $(date)
OUTPUT_FLAGSTAT=$FLAGSTAT_DIR/"$FILE"aligned_reads_sorted.flagstat
samtools flagstat $OUTPUT_SORTED_BAM -@ 24 -O tsv >$OUTPUT_FLAGSTAT

echo running samtools coverage $(date)
COVERAGE_DIR=$HOME_DIR/coverage
mkdir -p $COVERAGE_DIR
OUTPUT_COVERAGE=$COVERAGE_DIR/"$FILE"coverage.csv
OUTPUT_ABRIDGED_COVERAGE=$COVERAGE_DIR/"$FILE"abridged_coverage.csv
samtools coverage $OUTPUT_SORTED_BAM >$OUTPUT_COVERAGE
awk 'BEGIN{} {print $1","$6","$7 }' $OUTPUT_COVERAGE >$OUTPUT_ABRIDGED_COVERAGE

echo setting up for parallel processing of bcftools mpileup and call $(date)
CALL_DIR=$HOME_DIR/bcf_out
mkdir -p $CALL_DIR
TEMP_CALL_DIR=$CALL_DIR/"$FILE"temp
mkdir -p $TEMP_CALL_DIR
SEQUENCES=$(grep ">" $REF_INDEX | sed 's/>//' | cut -d" " -f1)

script_loc=/mnt/lz01/plachetzki/kcd1021/working_dir/913_project/bcftools_stuff.sh
#seq="$1"
#num=$2
#temp_call_dir=$3
#ref=$4
#infile=$5
set -f    # turn off globbing
IFS=$'\n' # split at newlines only
parallel -j 24 $script_loc {} '{#}' "$TEMP_CALL_DIR" "$REF_INDEX" "$OUTPUT_SORTED_BAM" ::: $SEQUENCES
unset IFS
set +f

shopt -s nullglob
BCF_FILES=($TEMP_CALL_DIR/*bcf)

echo concatenating bcf files together $(date)
bcftools concat ${BCF_FILES[@]} -o $CALL_DIR/"$FILE"complete.bcf &&
	rm -r $TEMP_CALL_DIR

echo done at $(date)
