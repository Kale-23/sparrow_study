#!/bin/bash

#SBATCH --job-name=TD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=/mnt/lz01/plachetzki/kcd1021/working_dir/913_project/output/TD/%j.out
#SBATCH --exclude=node\[117-118\]

##  Supply coverage directories space-separated in INDIRS , e.g. 
#INDIRS="premise_1x_cov/bcf_out/  premise_5x_cov/bcf_out/"  
#INDIRS="cov_1/bcf_out/ cov_3/bcf_out/ cov_5/bcf_out/ cov_7/bcf_out/ cov_10/bcf_out/ cov_15/bcf_out/ cov_20/bcf_out/"

## Alter ROOT for system-specific locality. 
#ROOT="~/working_dir/913_project/subs/" 

#usage() {
#	echo "You must specify the number of threads for bcftools merge"
#	exit 1
#}
#
#if [ $# -lt 1 ]; then
#	usage
#fi 
INDIRS=$1/bcf_out/
GREEDY=24 #threads
PWD=$(pwd) 

for dir in $INDIRS ; do  # Loop over all coverages
	FULLP=$dir
	echo $FULLP
	
	echo "0)"
	# 0) Get list of all files in this coverage directory
	BCFLIST=$FULLP"bcflist.txt"
	find $FULLP -type f -name "*.bcf" > $BCFLIST 
	sleep 1
	ALLFILES=$(cat $BCFLIST) 
	for elt in $ALLFILES; do
		echo $elt
	done

	echo "1)"
	# 1) All files must be compressed.
	cd $PWD
	cd $FULLP
	for basen in $ALLFILES; do
		if ! [ -f $basen".gz" ]; then
			bgzip $basen 
		fi
	done
	cd $PWD

	echo "2)"
	# 2) All files must be indexed
	cd $FULLP
	for basen in $ALLFILES; do
		if ! [ -f $basen".gz.csi" ]; then
			bcftools index --threads $GREEDY $basen".gz" 
		fi
	done
	cd $PWD 

	echo "3)"
	# 3) Merge all individuals into monolithic VCF
	BCFGZLIST=$FULLP"bcfgzlist.txt" 
	find $FULLP -type f -name "*.bcf.gz" > $BCFGZLIST # file-list for merge
	MONOVCF=$FULLP"mergemonolith.vcf.gz"
	if ! [ -f $MONOVCF ]; then
		bcftools merge --file-list $BCFGZLIST -Oz9 -o $MONOVCF --threads $GREEDY
	fi 

	echo "4)"
	#  4) Perform population statistic using the monolith.
	#    This is not multithreaded, but is luckily very fast. 
    OUTF=$FULLP"population_$(basename $(dirname $FULLP) | tr -d '/')"   # supplied output file.
	TAJIMA=$OUTF".Tajima.D"  # the output file actually created.
	vcftools --gzvcf $MONOVCF --TajimaD 71243 --out $OUTF

	echo "5)"
	#  5) Filter the nan's. Filter the contigs.
	cd $PWD
	cd $FULLP
	cat $TAJIMA | grep -v nan | grep -v NW_ > Dpopfilter

	echo "6)"
	#  6) Obtain a list of chromosome IDs. 
	CHROMOSOMES=$(cat Dpopfilter | grep NC_ | cut -f1 | sort | awk '!a[$0]++'  ) 

	echo "7)" 
	# 7) For each chromosome, create a .TD file that is just a list of Tajima D's one per line.
	for chr in $CHROMOSOMES; do
		OUTCHR=$FULLP$chr".TD"
		cat Dpopfilter | grep $chr | cut -f4 > $OUTCHR
	done
	cd $PWD

	# 8) Plot each of the .TD files with a box-whisker.  Should be approx 33 files.  
	#   Python or R. 
	#  You can skip steps 4 thru 7 since population_.Tajima.D is very simple format, 
	#    just make sure to filter NW_ contigs. 
done 

