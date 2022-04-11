#!/bin/bash

data_dir=$1 # full path to raw data dir (containing a separate folder with proper name for each sample)
adaptors=$2 # full path to fasta with adaptors/indexes
ref=$3 #dir with reference genome and plasmid (fasta, fa), no additional files

### MAIN ###

mkdir bowtie_ind
cd bowtie_ind
cat $ref* > reference
bowtie-build reference gen+p  # buiding common genome and plasmid bowtie indexes in advance
cd ../
mkdir gen_intervals
cd gen_intervals
python3 ../intervals.py --genome $ref*.fasta --interval 1000 # creating .bed file with genome intervals
cd ../

for dir in $(find $data_dir -mindepth 1 -type d ); do
    tag=$(basename $dir)
    raw="$(find $dir -type f -name '*.gz')" # .gz format only
    #bzcat $raw | gzip -c > ${tag}.gz   #can be used for re-compression from .bz2
    mkdir ${tag}_analysis
    cd ${tag}_analysis
    mkdir raw_QC
    fastqc -o raw_QC $raw
    echo "Trimmomatic started"
    trimmomatic SE $raw ${tag}_trimmed.fastq ILLUMINACLIP:$adaptors:2:30:10 SLIDINGWINDOW:3:20 CROP:24 MINLEN:14 # MAXLEN:21 can be also specified
    # cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 14 -M 24 -o ${tag}_trimmed.fastq *fastq.gz #alternative variant
    mkdir trimmed_QC
    fastqc -o trimmed_QC ${tag}_trimmed.fastq

    echo "bowtie alignment"
    bowtie -v 0 -m 1 -S --max ${tag}_multi.fq -x ../bowtie_ind/gen+p ${tag}_trimmed.fastq > ${tag}_align.sam # unique alignment
    bowtie -S -a --best --strata -v 0 -m 10000 -x ../bowtie_ind/gen+p ${tag}_multi.fq > ${tag}_multi.sam # multimappers
    
    echo "multimappers filtering"
    samtools view ${tag}_multi.sam | grep -P '\tNC_000913.3\t'| cut -f 1 | sort | uniq > multi_genome.txt # should be changed to BW25113 for yffP libraries
    samtools view ${tag}_multi.sam | grep -P '\tCloned_pBAD_CbuAg0_LacI_full\t' | cut -f 1 | sort | uniq > multi_plasmid.txt # should be changed for yffP libraries

    python3 ../filter_multimappers.py --mg multi_genome.txt --mp multi_plasmid.txt --multifq ${tag}_multi.fq # filtering and realignment
    bowtie -S -a --best --strata -v 0 -m 10000 -x ../bowtie_ind/gen+p uniq_multi.fq > ${tag}_u_multi.sam
    echo "SAMTOOLS"
    samtools view -S -b ${tag}_align.sam | samtools sort > ${tag}_align_sorted.bam
    samtools view -S -b ${tag}_u_multi.sam | samtools sort > ${tag}_u_multi_sorted.bam

    samtools merge -@ 8 ${tag}_final_sorted.bam *sorted.bam
    echo "alignment statistics"
    mkdir align_stat
    samtools view ${tag}_final_sorted.bam | grep -P  '\tNC_000913.3\t' | cut -f 1 | sort | uniq | wc -l > align_stat/aligned_on_g_p.txt   # should be changed to BW25113 for yffP libraries
    samtools view ${tag}_final_sorted.bam | grep -P  '\tCloned_pBAD_CbuAg0_LacI_full\t' | cut -f 1 | sort | uniq | wc -l >> align_stat/aligned_on_g_p.txt # should be changed for yffP libraries

    python3 ../alignment_stat.py --genome $ref*.fasta --plasmid $ref*.fa --align align_stat/aligned_on_g_p.txt --copy_num 12 # pBAD copy number from literature
    mv alignment_stat.txt align_stat/alignment_stat.txt

    echo "intersecting genome intervals, counting coverage"
    bedtools intersect -a ../gen_intervals/intervals*.bed -b ${tag}_final_sorted.bam -wa -wb -F 0.5 > ${tag}_intersected.tsv
    cut -f 8 ${tag}_intersected.tsv | sort | uniq -c > ${tag}_counts.txt
    python3 ../count_coverage.py --intersect ${tag}_intersected.tsv --counts ${tag}_counts.txt

    samtools view -F 16 -b ${tag}_final_sorted.bam > ${tag}_plus.bam # repeating previous steps for + strand
    bedtools intersect -a ../gen_intervals/intervals*.bed -b ${tag}_plus.bam -wa -wb -F 0.5 > ${tag}plus_intersected.tsv
    python3 ../count_coverage.py --intersect ${tag}plus_intersected.tsv --counts ${tag}_counts.txt

    samtools view -f 16 -b ${tag}_final_sorted.bam > ${tag}_minus.bam # repeating previous steps for - strand
    bedtools intersect -a ../gen_intervals/intervals*.bed -b ${tag}_minus.bam -wa -wb -F 0.5 > ${tag}minus_intersected.tsv
    python3 ../count_coverage.py --intersect ${tag}minus_intersected.tsv --counts ${tag}_counts.txt
    echo "COMPLETED"
    cd ..
done




