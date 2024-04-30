#!/bin/bash

in_bam=$1
out_dir=$2
out_forward=$out_dir/$(basename $in_bam)_R1.fastq.gz
out_reverse=$out_dir/$(basename $in_bam)_R2.fastq.gz
in_file=${in_bam}_aligned_reads_sorted.bam
samtools collate -@10 -u -O $in_file | \
samtools fastq -@10 -c6 -1 $out_forward -2 $out_reverse -0 /dev/null -s /dev/null -n

output_dir=${out_dir%%mito_reads/}de_novo_alignment/$(basename $in_bam)
mkdir -p $output_dir

spades.py \
    -1 $out_forward \
    -2 $out_reverse \
    -o $output_dir \
    -m 250 \
    --plasmid \
    -t 24
