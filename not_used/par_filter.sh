#!/usr/bin/env bash

database=$1
query=$2
nodes=$3
user=$4

usage="
USAGE:\n\tlaunch_blast <database> <query> <nodes> <user>\n
\t<database>:\tBLAST database version 5\n
\t<query>:\t*.list file obtained from par_split\n
"

# I check if user is defined
[ -z $database  ] && echo "Database name is missing" && echo -e $usage
[ -z $database  ] && exit 1
[ -z $query  ] && echo "Query is missing!" && echo -e $usage && exit 1
[ -z $nodes  ] && echo "Which nodes do you want to use?" && echo -e $usage && exit 1
[ -z $user  ] && echo "Username is missing" && echo -e $usage
[ -z $user  ] && exit 1

#echo "Check point ...[OK]"
#exit 0
# Check point ...[OK]

qpath=$(dirname $query)
qname=$(basename $query)
[ $qpath == "."  ] && silpath="sil.txt"
[ $qpath != "."  ] && [ $qpath != ""  ] && silpath="${qpath}/sil.txt"

#echo $dbpath
#echo $dbname
#echo $silpath
#echo "Check point ...[OK]"
#exit 0
# Check point successes
# Check point ...[OK]

#	Create sil.txt in shared_db folder.
#	From now on, you need to test it on server master-006.

load_module ncbi-blast/2.13.0+/bin
blastdbcmd -db ${shared_db}/${database} -entry all -outfmt %a -out $silpath
#	Parallel execution.

#echo "Check point ...[OK]"
# Now, sil.txt file is inside the data path. And it is needed to transfer the file to node.
# So, how can I transfer this file in every node?
# Let me think!
#exit 0

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
		-task blastn \
		-out {1.}.blast \
		-outfmt 6 \
		-evalue 50 \
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
		-x ${shared_db}/random/ecoli \
		-fU {1.}.outseq.fasta \
		-S {1.}.outseq.sam; \
	samtools view -f 0 -b {1.}.outseq.sam > {1.}.outseq.bam; \
	samtools sort -l 9 -O BAM {1.}.outseq.bam > {1.}.outseq.sort.bam; \
	samtools index -b {1.}.outseq.sort.bam > {1.}.outseq.sort.bam.bai; \
	stringtie {1.}.outseq.sort.bam -o {1.}.outseq.gtf -L -s 1 -g 250" ::: $silpath
