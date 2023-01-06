#!/bin/bash

### Download MSU7.0 cDNA
#wget ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.cdna


## Map reads to CDS using Bowtie2
bowtie2-build --offrate 1 all.cdna MSU7_cDNA


### Rep1
for i in APO IR64 F1; do
( bowtie2 --sensitive --phred33 -a --time --threads 24 --maxins 350 -x MSU7_cDNA -1 Rep1_${i}-C_1.fq -2 Rep1_${i}-C_2.fq | samtools view -bS - > Rep1_${i}-C.bam ) > Rep1_${i}-C_bowtie2.log 2>&1
( bowtie2 --sensitive --phred33 -a --time --threads 24 --maxins 350 -x MSU7_cDNA -1 Rep1_${i}-S_1.fq -2 Rep1_${i}-S_2.fq | samtools view -bS - > Rep1_${i}-S.bam ) > Rep1_${i}-S_bowtie2.log 2>&1
done


### Rep2
for i in APO IR64 F1; do
( bowtie2 --sensitive --phred64 -a --time --threads 24 --maxins 500 -x MSU7_cDNA -1 Rep2_${i}-C_1.fq -2 Rep2_${i}-C_2.fq | samtools view -bS - > Rep2_${i}-C.bam ) > Rep2_${i}-C_bowtie2.log 2>&1
( bowtie2 --sensitive --phred64 -a --time --threads 24 --maxins 500 -x MSU7_cDNA -1 Rep2_${i}-S_1.fq -2 Rep2_${i}-S_2.fq | samtools view -bS - > Rep2_${i}-S.bam ) > Rep2_${i}-S_bowtie2.log 2>&1
done



