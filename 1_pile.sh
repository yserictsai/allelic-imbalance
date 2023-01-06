#!/bin/bash

samtools mpileup -I -f MSU7_all.cdna Rep1_APO-C_sorted.bam Rep1_APO-S_sorted.bam Rep1_IR64-C_sorted.bam Rep1_IR64-S_sorted.bam  > Rep1_pool_mpileup.txt &
samtools mpileup -I -f MSU7_all.cdna Rep2_APO-C_sorted.bam Rep2_APO-S_sorted.bam Rep2_IR64-C_sorted.bam Rep2_IR64-S_sorted.bam  > Rep2_pool_mpileup.txt &
