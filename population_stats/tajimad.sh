#!/bin/bash
##  Supply coverage directories space-separated in INDIRS , e.g. 
#INDIRS="premise_1x_cov/bcf_out/  premise_5x_cov/bcf_out/"  
INDIRS="premise_1x_cov/bcf_out/sandbox/"

# Alter ROOT for system-specific locality. 
ROOT="/home/unhTW/share/mcbs913_2024/group3/" 

usage() {
	echo "You must specify the number of threads for bcftools merge"
	exit 1
}

if [ $# -lt 1 ]; then
	usage
fi 
GREEDY=$1
PWD=$(pwd) 

for dir in $INDIRS ; do  # Loop over all coverages
	FULLP=$ROOT$dir
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
			bcftools index $basen".gz" 
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
	OUTF=$FULLP"population_"   # supplied output file.
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

