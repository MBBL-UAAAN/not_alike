#!/usr/bin/env bash

OPTIND=1
usage="
USAGE:\n
	not_alike [arguments]\n
			-i|--input <Query fasta file>\n
			-np|--num-procs <Number of CPU cores>\n
			-n|--nodes <Nodes names>\n
			-fs|--fragment-size <Fragment size in nt>\n
			-ss|--step-size <Step size in nt>\n
			-b|--batch-number <Number of batch query files>\n
			-d|--database <Full path name of database file>\n
			-hi|--hisat2-idx <Full path name of hisat2 index files>\n
			-id|--identity <Query percentage of identity cutoff>\n
			-qc|--qcov <Query HSP coverage percenatge cutoff>\n
			-t|--task <Blastn task [blastn, megablast, dc-megablast]>\n
			-e|--evaue <E-value cutoff>
"

#	not_alike-pipeline general arguments
input=""
nproc=""
nodes=""
frg_size=""
stp_size=""
batch_num=""
database=""

#	Blastn specific arguments
task=""
ident=""
qcov=""
evalue=""

while [[ "$#" -gt 0  ]];
do

	case $1 in
		-i|--input) input=$2; shift ;;
		-np|--num-procs) nproc=$2; shift ;;
		-n|--nodes) nodes=$2; shift ;;
		-fs|--fragment-size) frg_size=$2; shift ;;
		-ss|--step-size) stp_size=$2; shift ;;
		-b|--batch-number) batch_num=$2; shift ;;
		-d|--database) database=$2; shift ;;
		-hi|--hisat2-idx) ht2_idx=$2; shift ;;
		-id|--identity) ident=$2; shift ;;
		-qc|--qcov) qcov=$2; shift ;;
		-t|--task) task=$2; shift ;;
		-e|--evalue) evalue=$2; shift ;;
		*) echo "Unknown parameter passed $1"; echo -e $usage; exit 1 ;;
	esac
	shift
done

[ -z $input  ] && echo "Input query fasta file is empty" && echo -e $usage && exit 1
[ -z $nproc  ] && echo "Number of processors is empty" && echo -e $usage && exit 1
[ -z $nodes  ] && echo "No nodes were defined" && echo -e $usage && exit 1
[ -z $frg_size  ] && echo "Fragment size argument is empty" && echo -e $usage && exit 1
[ -z $stp_size  ] && echo "Step size argument is empty" && echo -e $usage && exit 1
[ -z $batch_num  ] && echo "Number of batch files argument is empty" && echo -e $usage && exit 1
[ -z $database ] && echo "Database full path name is empty" && echo -e $usage && exit 1
[ -z $ht2_idx ] && echo "Hisat2 index file full path name is empty" && echo -e $usage && exit 1
[ -z $ident ] && echo "Percentage of identity was not defined" && echo -e $usage && exit 1
[ -z $qcov ] && echo "Query HSP coverage percentage was not defined" && echo -e $usage && exit 1
[ -z $task ] && echo "Blastn task was not defined" && echo -e $usage && exit 1
[ -z $evalue ] && echo "E-value was not defined" && echo -e $usage && exit 1

infile=$(basename $input)
route=$(dirname $input)
PID=${RANDOM}${RANDOM}

parallel -j $nproc \
	--env load_module,shared \
	--sshlogin $nodes \
	--transferfile {} \
	--cleanup \
	--return {.}_split.tar.gz \
	--cleanup \
	"load_module parallel_useq/0.0.1/scripts; \
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

query=$route/${infile%.*}.list
qpath=$(dirname $query)
qname=$(basename $query)
[ $qpath == "."  ] && silpath="sil.txt"
[ $qpath != "."  ] && [ $qpath != ""  ] && silpath="${qpath}/sil.txt"

load_module ncbi-blast/2.13.0+/bin
blastdbcmd -db ${shared_db}/${database} -entry all -outfmt %a -out $silpath

parallel -j 1 \
	--env load_module,shared,shared_db \
	--sshlogin $nodes \
	--transferfile {1} \
	--transferfile {2} \
	--return {1.}.outseq.gtf \
	--cleanup \
	-a $query \
	"load_module parallel_useq/0.0.1/scripts; \
	load_module ncbi-blast/2.13.0+/bin; \
	load_module hisat2/2.2.1; \
	load_module samtools/1.15.1; \
	load_module stringtie/2.2.1; \
	echo \"Working on {1}\"; \
	grep '>' {1} | tr -d '/^>/' > {1.}.guid.txt; \
	cat {2} | \
	while read -r line; \
	do
		echo $line > {1//}/sid.txt; \
		blastn -db ${shared_db}/${database} \
		-query {1} \
		-seqidlist {1//}/sid.txt \
		-task $task \
		-max_target_seqs 5 \
		-perc_identity $ident \
		-qcov_hsp_perc $qcov \
		-out {1.}.blast \
		-outfmt 6 \
		-evalue $evalue \
		-gapopen 0 \
		-gapextend 2 \
		-penalty -1 \
		-reward 1; \
		cut -f1 {1.}.blast | sort | uniq  > {1.}.guid_hit.txt; \
		update-seqs {1} {1.}.guid.txt {1.}.guid_hit.txt > {1.}.tmp; \
		mv {1.}.tmp {1}; \
		[ ! -s {1}  ] && break; \
	done; \
	mv {1} {1.}.outseq.fasta; \
	hisat2 --no-spliced-alignment \
		--no-templatelen-adjustment \
		-x ${shared_db}/$ht2_idx \
		-fU {1.}.outseq.fasta \
		-S {1.}.outseq.sam; \
	samtools view -f 0 -b {1.}.outseq.sam > {1.}.outseq.bam; \
	samtools sort -l 9 -O BAM {1.}.outseq.bam > {1.}.outseq.sort.bam; \
	samtools index -b {1.}.outseq.sort.bam > {1.}.outseq.sort.bam.bai; \
	stringtie {1.}.outseq.sort.bam -o {1.}.outseq.gtf -L -s 1 -g 250" ::: $silpath
