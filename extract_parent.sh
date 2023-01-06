#!/bin/bash

perl extract_apo_1c.pl mod_Rep1_APO-C.bam all_rep_snpcall.txt &
perl extract_apo_1s.pl mod_Rep1_APO-S.bam all_rep_snpcall.txt &
perl extract_apo_2c.pl mod_Rep2_APO-C.bam all_rep_snpcall.txt &
perl extract_apo_2s.pl mod_Rep2_APO-S.bam all_rep_snpcall.txt &

perl extract_ir64_1c.pl mod_Rep1_IR64-C.bam all_rep_snpcall.txt &
perl extract_ir64_1s.pl mod_Rep1_IR64-S.bam all_rep_snpcall.txt &
perl extract_ir64_2c.pl mod_Rep2_IR64-C.bam all_rep_snpcall.txt &
perl extract_ir64_2s.pl mod_Rep2_IR64-S.bam all_rep_snpcall.txt &




