#!/bin/bash

for i in APO IR64;do
samtools sort mod_Rep2_${i}-S.bam mod_Rep2_${i}-S_sorted &
samtools sort mod_Rep2_${i}-C.bam mod_Rep2_${i}-C_sorted &
done




for i in APO IR64;do
samtools sort mod_Rep1_${i}-S.bam mod_Rep1_${i}-S_sorted &
samtools sort mod_Rep1_${i}-C.bam mod_Rep1_${i}-C_sorted &
done
