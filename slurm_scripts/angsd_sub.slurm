#! /bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --job-name="ANGSD"
#SBATCH --mem 150G
#SBATCH --output=ANGSD_%j.log
#SBATCH --exclude=node117,node118

bam_file=$1
new_dir=$2
angsd -GL 1 \
    -out $new_dir/genolike \
    -nThreads 20 \
    -doGlf 2 \
    -doMajorMinor 1  \
    -doMaf 2 \
    -SNP_pval 1e-6 \
    -bam $bam_file 
