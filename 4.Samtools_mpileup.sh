#!/bin/bash


for i in `ls *.bam | sed 's/.bam//'`; do
samtools sort ${i}.bam ${i}_sorted
done


for i in Rep1 Rep2; do
samtools mpileup -I -f all.cds ${i}_APO-C_sorted.bam ${i}_IR64-C_sorted.bam ${i}_F1-C_sorted.bam | ./Parse_mpileup.pl > Parsed_${i}-C_combined.mpileup
samtools mpileup -I -f all.cds ${i}_APO-S_sorted.bam ${i}_IR64-S_sorted.bam ${i}_F1-S_sorted.bam | ./Parse_mpileup.pl > Parsed_${i}-S_combined.mpileup
done
