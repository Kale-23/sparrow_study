#! /bin/bash

shopt -s nullglob
for dir in /mnt/lz01/plachetzki/kcd1021/working_dir/913_project/mito/*; do
    dir_name=$(basename $dir)
    new_dir=$dir/pcangsd
    mkdir -p $new_dir
    file=$dir/angsd/genolike.beagle.gz
    output=$new_dir/pcangsd_${dir_name}_output
    sbatch ./pcangsd.slurm $file $output 
done
