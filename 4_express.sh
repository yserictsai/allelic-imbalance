#!/usr/bin/bash


for i in APO IR64 F1; do 
express --output-dir express_Rep1_${i}_C modwithindel_ref.fasta mod_Rep1_${i}-C.bam & 
express --output-dir express_Rep1_${i}_S modwithindel_ref.fasta mod_Rep1_${i}-S.bam & 
done


for i in APO IR64 F1; do
express --output-dir express_Rep2_${i}_C modwithindel_ref.fasta mod_Rep2_${i}-C.bam &
express --output-dir express_Rep2_${i}_S modwithindel_ref.fasta mod_Rep2_${i}-S.bam &
done
