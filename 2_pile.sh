#!/bin/bash

samtools mpileup -I -f modwithindel_ref.fasta mod_Rep1_APO-C_sorted.bam mod_Rep1_APO-S_sorted.bam mod_Rep1_IR64-C_sorted.bam mod_Rep1_IR64-S_sorted.bam  > mod_Rep1_pool_mpileup.txt &
samtools mpileup -I -f modwithindel_ref.fasta mod_Rep2_APO-C_sorted.bam mod_Rep2_APO-S_sorted.bam mod_Rep2_IR64-C_sorted.bam mod_Rep2_IR64-S_sorted.bam  > mod_Rep2_pool_mpileup.txt &
