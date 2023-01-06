#!/bin/bash

for i in `ls Rep1*.fq | sed s/.fq//g`; do
fastq_quality_filter -i ${i}.fq -q30 -p80 -Q33 | fastq_quality_trimmer -t20 -l15 -Q33 -v > QC_${i}.fq &
done

for j in `ls Rep2*.fq | sed s/.fq//g`; do
fastq_quality_filter -i ${j}.fq -q30 -p80 | fastq_quality_trimmer -t20 -l15 -v > QC_${j}.fq &
done
