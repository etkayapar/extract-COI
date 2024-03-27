#! /usr/bin/env bash

# ------------------------------------------------------------------
# Usage:
#       ./main.sh <path/to/genome/paths/file.txt> <NUM_THREADS>
# ------------------------------------------------------------------

genome_paths=$1 ## path to a text file with full paths to assemblies
threads=$2 ## number of independent parallel processes to run

scr_root=$(dirname "$0")
scr_root=$(realpath "$scr_root")

## Create blast databases ------------------------------------------
mkdir -p databases
if [[ ! -n $(ls -A databases/) ]]
then
	rm -rf blast_results/ all_COI.fa*
	echo -e "Creating blast databases for the input assemblies..."
	cat ${genome_paths} | xargs -n 1 -P $threads ${scr_root}/makeblastdb_single.sh
	echo -e "Done.\n"
else
	echo -e "Skipping blast db creation since './databases/' is not empty\n"
fi

## Run blast -------------------------------------------------------
mkdir -p blast_results
if [[ ! -n $(ls -A blast_results | grep "blastn.tsv") ]]
then
	rm -rf blast_results/*.fasta all_COI.fa*
	echo -e "Running blast..."
	cat ${genome_paths} | xargs -n1 -P $threads ${scr_root}/run_blast_single.sh
	echo -e "Done.\n"
else
	echo -e "Skipping running blast since './blast_results/' has results from previous runs\n"
fi


## Extract best COI hits -------------------------------------------
if [[ ! -n $(ls -A blast_results | grep -E "_COI_.+\.fasta") ]]
then
	rm -rf all_COI.fa*
	echo -e "Extracting COI sequences according to blast results..."
	cp ${genome_paths} blast_results/genomes.txt
	cd blast_results
	${scr_root}/extract_seq_from_blast.py
	cd ..
	echo -e "Done.\n"
else
	echo -e "Skipping running the COI sequence extraction step.\n"
fi

cat blast_results/*_COI_nt.fasta > all_COI.fa
cat blast_results/*_COI_aa.fasta > all_COI.faa
