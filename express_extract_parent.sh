#!/usr/bin/bash


#samtools view -t modwithindel_ref.fasta.fai -Sb mod_Rep1_extract_APO_c.sam | express --output-dir rep1_express_e_apo_c modwithindel_ref.fasta - &
samtools view -t modwithindel_ref.fasta.fai -Sb mod_Rep2_extract_APO_c.sam | express --output-dir rep2_express_e_apo_c modwithindel_ref.fasta - &
samtools view -t modwithindel_ref.fasta.fai -Sb mod_Rep1_extract_APO_s.sam | express --output-dir rep1_express_e_apo_s modwithindel_ref.fasta - &
samtools view -t modwithindel_ref.fasta.fai -Sb mod_Rep2_extract_APO_s.sam | express --output-dir rep2_express_e_apo_s modwithindel_ref.fasta - &
samtools view -t modwithindel_ref.fasta.fai -Sb mod_Rep1_extract_IR64_c.sam | express --output-dir rep1_express_e_ir64_c modwithindel_ref.fasta - &
samtools view -t modwithindel_ref.fasta.fai -Sb mod_Rep2_extract_IR64_c.sam | express --output-dir rep2_express_e_ir64_c modwithindel_ref.fasta - &
samtools view -t modwithindel_ref.fasta.fai -Sb mod_Rep1_extract_IR64_s.sam | express --output-dir rep1_express_e_ir64_s modwithindel_ref.fasta - &
samtools view -t modwithindel_ref.fasta.fai -Sb mod_Rep1_extract_IR64_s.sam | express --output-dir rep2_express_e_ir64_s modwithindel_ref.fasta - &



