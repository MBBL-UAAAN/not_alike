#!/usr/bin/env bash

###	Parallel split genome.


input=$1
nproc=$2
nodes=$3
frg_size=$4
stp_size=$5
batch_num=$6

[ -z $input  ] && echo "Input is a first place argument, but it is missing!"
[ -z $nproc  ] && echo "Number of processors is a second place argument, but it is missing!"
[ -z $nodes  ] && echo "Number of nodes is a third place argument, but it is missing!"
[ -z $frg_size  ] && echo "Fragment size is a fourth place argument, but it is missing!"
[ -z $stp_size  ] && echo "Step size is a fifth place argument, but it is missing!"
[ -z $batch_num  ] && echo "Batch number is a sixth place argument, but it is missing!"
[ -z $input  ] && exit 1
[ -z $nproc  ] && exit 1
[ -z $nodes  ] && exit 1
[ -z $frg_size  ] && exit 1
[ -z $stp_size  ] && exit 1
[ -z $batch_num  ] && exit 1

infile=$(basename $input)
route=$(dirname $input)
PID=${RANDOM}${RANDOM}

parallel -j $nproc \
	--env load_module,shared \
	--sshlogin $nodes \
	--joblog split_genome_$PID \
	--transferfile {} \
	--cleanup \
	--return {.}_split.tar.gz \
	--cleanup \
	"load_module parallel_blast/scripts; \
	echo \"Sending split-genome process to $nodes\";\
	split-genome {} $frg_size $stp_size {.}_split 1 $batch_num; \
	tar --gzip -cvf {.}_split.tar.gz {.}_split*" ::: $input

echo "Extracting ${input%.*}_split.tar.gz"
tar -zxvf ${input%.*}_split.tar.gz
echo "Creating a list file with input file names to the next prosses"
for file in $route/*_split_*;
do
	echo $file >> $route/${infile%.*}.list
done


