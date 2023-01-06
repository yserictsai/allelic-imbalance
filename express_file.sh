#!/usr/bin/bash


for i in APO IR64 F1; do

less express_Rep1_${i}_C/results.xprs | cut -f2,5 > express_Rep1_${i}_C_tot_count.txt &
less express_Rep1_${i}_S/results.xprs | cut -f2,5 > express_Rep1_${i}_S_tot_count.txt &
less express_Rep2_${i}_C/results.xprs | cut -f2,5 > express_Rep2_${i}_C_tot_count.txt &
less express_Rep2_${i}_S/results.xprs | cut -f2,5 > express_Rep2_${i}_S_tot_count.txt &



done



for i in apo ir64; do

less rep1_express_f1_from_${i}_c/results.xprs | cut -f2,5 > express_Rep1_from_${i}_C_tot_count.txt &
less rep1_express_f1_from_${i}_s/results.xprs | cut -f2,5 > express_Rep1_from_${i}_S_tot_count.txt &
less rep2_express_f1_from_${i}_c/results.xprs | cut -f2,5 > express_Rep2_from_${i}_C_tot_count.txt &
less rep2_express_f1_from_${i}_s/results.xprs | cut -f2,5 > express_Rep2_from_${i}_S_tot_count.txt &


done



for i in rep1 rep2; do

less ${i}_express_e_apo_c/results.xprs | cut -f2,5 > ${i}_e_apo_c_tot_count.txt &
less ${i}_express_e_apo_s/results.xprs | cut -f2,5 > ${i}_e_apo_s_tot_count.txt &
less ${i}_express_e_ir64_c/results.xprs | cut -f2,5 > ${i}_e_ir64_c_tot_count.txt &
less ${i}_express_e_ir64_s/results.xprs | cut -f2,5 > ${i}_e_ir64_s_tot_count.txt &


done
