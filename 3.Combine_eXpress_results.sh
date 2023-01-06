#!/bin/bash


### Rep1
for i in APO IR64 F1; do
grep -v '^bundle_id' Rep1_${i}-C_eXpress/results.xprs | sort -k2,7 | cut -f2,8,11 > Rep1_${i}-C_eXpress.temp
grep -v '^bundle_id' Rep1_${i}-S_eXpress/results.xprs | sort -k2,7 | cut -f2,8,11 > Rep1_${i}-S_eXpress.temp
done
echo "CDS\tAPO-C_counts\tAPO-C_FPKM\tIR64-C_counts\tIR64-C_FPKM\tF1-C_counts\tF1-C_FPKM" > Rep1-C_counts_FPKM.tab
paste Rep1_APO-C_eXpress.temp Rep1_IR64-C_eXpress.temp Rep1_F1-C_eXpress.temp | cut -f1,2,3,5,6,8,9 >> Rep1-C_counts_FPKM.tab
echo "CDS\tAPO-S_counts\tAPO-S_FPKM\tIR64-S_counts\tIR64-S_FPKM\tF1-S_counts\tF1-S_FPKM" > Rep1-S_counts_FPKM.tab
paste Rep1_APO-S_eXpress.temp Rep1_IR64-S_eXpress.temp Rep1_F1-S_eXpress.temp | cut -f1,2,3,5,6,8,9 >> Rep1-S_counts_FPKM.tab
rm -f *.temp


### Rep2
for i in APO IR64 F1; do
grep -v '^bundle_id' Rep2_${i}-C_eXpress/results.xprs | sort -k2,7 | cut -f2,8,11 > Rep2_${i}-C_eXpress.temp
grep -v '^bundle_id' Rep2_${i}-S_eXpress/results.xprs | sort -k2,7 | cut -f2,8,11 > Rep2_${i}-S_eXpress.temp
done
echo "CDS\tAPO-C_counts\tAPO-C_FPKM\tIR64-C_counts\tIR64-C_FPKM\tF1-C_counts\tF1-C_FPKM" > Rep2-C_counts_FPKM.tab
paste Rep2_APO-C_eXpress.temp Rep2_IR64-C_eXpress.temp Rep2_F1-C_eXpress.temp | cut -f1,2,3,5,6,8,9 >> Rep2-C_counts_FPKM.tab
echo "CDS\tAPO-S_counts\tAPO-S_FPKM\tIR64-S_counts\tIR64-S_FPKM\tF1-S_counts\tF1-S_FPKM" > Rep2-S_counts_FPKM.tab
paste Rep2_APO-S_eXpress.temp Rep2_IR64-S_eXpress.temp Rep2_F1-S_eXpress.temp | cut -f1,2,3,5,6,8,9 >> Rep2-S_counts_FPKM.tab
rm -f *.temp

