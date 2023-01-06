#!/bin/bash


for i in Rep1 Rep2; do
perl Combine_FPKM_SNPcounts.pl ${i}-C_counts_FPKM.tab Parsed_${i}-C_combined.mpileup > Final_${i}-C.csv
perl Combine_FPKM_SNPcounts.pl ${i}-S_counts_FPKM.tab Parsed_${i}-S_combined.mpileup > Final_${i}-S.csv
done
