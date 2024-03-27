#!/bin/bash

sp=$1
sp=${sp##*/}
sp=${sp%%.*}

tblastn -db databases/${sp} -query COI_query.fasta -out blast_results/${sp}_tblastn.tsv -task tblastn -outfmt "7 sseqid qseqid sstart send sframe evalue bitscore qlen length"
