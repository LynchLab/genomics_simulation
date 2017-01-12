#!/bin/bash
reference=$1
name=$2
out=$3
ID=$4

bwa index $1

bwa aln $reference $name.1.fq > $name.1.sai
bwa aln $reference $name.2.fq > $name.2.sai

bwa sampe $reference $name.1.sai $name.2.sai $name.1.fq $name.2.fq -r "@RG\tID:$ID\tSM:$ID" | samtools view -S -b - | samtools sort - | samtools rmdup - - > $out.sort.rmdup.bam

samtools faidx $reference 
samtools view -H $out.sort.rmdup.bam > $name-header.txt
