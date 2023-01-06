#!/bin/bash


### Rep1
for i in APO IR64 F1; do
express --output-dir Rep1_${i}-C_eXpress all.cds Rep1_${i}-C.bam &
express --output-dir Rep1_${i}-S_eXpress all.cds Rep1_${i}-S.bam &
done


### Rep2
for i in APO IR64 F1; do
express --output-dir Rep2_${i}-C_eXpress all.cds Rep2_${i}-C.bam &
express --output-dir Rep2_${i}-S_eXpress all.cds Rep2_${i}-S.bam &
done
