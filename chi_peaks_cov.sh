#!/bin/bash

analysis=$1 #path to dir with previous analysis
chi_bed=$2 #path to gene-chi intervals bed (e.g second left chi - gene start, gene end - second right chi)
tag=$3 #sample name

### MAIN ###

echo 'Computing coverage of the regions between the second Chi-sites'
mkdir ${tag}_chi_peaks
cd ${tag}_chi_peaks

bedtools intersect -a $chi_bed -b $analysis/*final_sorted.bam -wa -wb -F 0.51 > 1_table_intersected.tsv

samtools view $analysis/*final_sorted.bam | grep -P "\tNC_" | cut -f 1 | sort | uniq -c > counts.txt # or \tBW
python3 ../count_coverage.py --intersect 1_table_intersected.tsv --counts counts.txt

cd ..
echo 'DONE'

