#! /bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name="pcangsd"
#SBATCH --mem 150G
#SBATCH --output=pcangsd_%j.log
#SBATCH --exclude=node117,node118

input=$1
output=$2

pcangsd -b $input \
    -o $output \
    -e 15 \
    --selection \
    --admix \
    --tree \
    --maf_save \
    --pi_save \
    --dosage_save \
    --sites_save \
    -t 24
