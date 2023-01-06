#!/bin/bash


for i in APO IR64 F1;do
bowtie2 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4  -X 500 --no-mixed --no-discordant --fr -x mod_indel_MSU7 -q --phred33 -p 20 -1 ../FastQC/Rep1_${i}-C_1.fq  -2 ../FastQC/Rep1_${i}-C_2.fq | samtools view -bS - -o mod_Rep1_${i}-C.bam 
bowtie2 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4  -X 500 --no-mixed --no-discordant --fr -x mod_indel_MSU7 -q --phred33 -p 20 -1 ../FastQC/Rep1_${i}-S_1.fq  -2 ../FastQC/Rep1_${i}-S_2.fq | samtools view -bS - -o mod_Rep1_${i}-S.bam 
done



for i in APO IR64 F1;do
bowtie2 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4  -X 600 --no-mixed --no-discordant --fr -x mod_indel_MSU7 -q --phred64 -p 20 -1 ../FastQC/QC_Rep2_${i}-C_1.fq  -2 ../FastQC/QC_Rep2_${i}-C_2.fq | samtools view -bS - -o mod_Rep2_${i}-C.bam 
bowtie2 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4  -X 600 --no-mixed --no-discordant --fr -x mod_indel_MSU7 -q --phred64 -p 20 -1 ../FastQC/QC_Rep2_${i}-S_1.fq  -2 ../FastQC/QC_Rep2_${i}-S_2.fq | samtools view -bS - -o mod_Rep2_${i}-S.bam 
done
