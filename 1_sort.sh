
#!/bin/bash

for i in APO IR64 F1;do
samtools sort Rep2_${i}-S.bam Rep2_${i}-S_sorted &
samtools sort Rep2_${i}-C.bam Rep2_${i}-C_sorted &
done




for i in APO IR64 F1;do
samtools sort Rep1_${i}-S.bam Rep1_${i}-S_sorted &
samtools sort Rep1_${i}-C.bam Rep1_${i}-C_sorted &
done

