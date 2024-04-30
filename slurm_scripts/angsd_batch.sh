#! /bin/bash
shopt -s nullglob
for dir in /mnt/lz01/plachetzki/kcd1021/working_dir/913_project/mito/*; do
    dir_name=$(basename $dir)
    new_dir=$dir/angsd
    mkdir -p $new_dir
    bam_list=$new_dir/bam_list_${dir_name}.txt
    rm $bam_list
    for item in $dir/alignment_out/*.bam; do 
        echo $item >> $bam_list
    done
    sbatch ./angsd_sub.slurm $bam_list $new_dir
done
