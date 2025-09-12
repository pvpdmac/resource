#!/bin/bash
set -e

dir_work=/path/to/data
dir_fastq=${dir_work}/Fastq
dir_genome=/path/to/genome

samples=${dir_work}/samples.txt
gtf_file=${dir_genome}/genome.gtf
fasta_file=${dir_genome}/genome.fna


# FastQC
# This will run quietly in the background by default; remove the trailing "&" to run it in the foreground, sequentially
if [ ! -d FastQC ]; then mkdir FastQC; fi
fastqc -t 6 -q -o FastQC/ ${dir_fastq}/*.fastq.gz &
p1=$!


# STAR
if [ ! -d STAR ]; then mkdir STAR; fi
while read sample; do
  STAR \
    --runMode alignReads \
    --runThreadN 16 \
    --genomeDir ${dir_genome} \
    --genomeLoad LoadAndKeep \
    --limitBAMsortRAM 48000000000 \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesIn ${dir_fastq}/${sample}_R1.fastq.gz ${dir_fastq}/${sample}_R2.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix STAR/${sample}_
done < ${samples}

rename 's/_Aligned.sortedByCoord.out//' STAR/*.bam

STAR --genomeDir ${dir_genome} --genomeLoad Remove
rm Aligned.out.sam Log.out Log.final.out Log.progress.out SJ.out.tab


# HTSeq
if [ ! -d HTSeq ]; then mkdir HTSeq; fi
find STAR/ -name "*.bam" | parallel --jobs 6 \
  "htseq-count -s reverse -a 10 -r pos {} "${gtf_file}" > HTSeq/{/.}.count"


# MultiQC
# Wait for FastQC to finish before running MultiQC
wait $p1
multiqc \
  -f \
  -o MultiQC \
  FastQC/ \
  STAR/ \
  HTSeq/
