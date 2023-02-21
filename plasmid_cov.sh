#!/bin/bash

tag=$1 # lib name
ref=$2 # full path to fasta with reference plasmid

### MAIN ###



bowtie-build $ref ind


bowtie -S -x ind ../../new_ago_2022/${tag}_analysis/${tag}_trimmed.fastq > ${tag}_align.sam
samtools view -S -b ${tag}_align.sam | samtools sort > ${tag}_align_sorted.bam
samtools view -F 16 -b ${tag}_align_sorted.bam > ${tag}_plus.bam
samtools view -f 16 -b ${tag}_align_sorted.bam > ${tag}_minus.bam
samtools depth ${tag}_align_sorted.bam > ${tag}_p.coverage
samtools depth ${tag}_plus.bam > ${tag}_p_plus.coverage
samtools depth ${tag}_minus.bam > ${tag}_p_minus.coverage

echo "completed"
