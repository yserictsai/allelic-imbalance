#!/bin/bash

samtools view -h mod_Rep1_F1-C.bam -o mod_Rep1_F1-C.sam
samtools view -h mod_Rep1_F1-S.bam -o mod_Rep1_F1-S.sam
samtools view -h mod_Rep2_F1-C.bam -o mod_Rep1_F1-C.sam
samtools view -h mod_Rep2_F1-S.bam -o mod_Rep2_F1-S.sam

perl extract_f1.pl all_rep_snpcall.txt mod_Rep1_F1-C.sam 
perl extract_f1.pl all_rep_snpcall.txt mod_Rep1_F1-S.sam
perl extract_f1.pl all_rep_snpcall.txt mod_Rep1_F1-C.sam
perl extract_f1.pl all_rep_snpcall.txt mod_Rep2_F1-S.sam
