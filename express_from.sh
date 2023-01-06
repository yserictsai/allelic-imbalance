#!/usr/bin/bash


samtools view -t modwithindel_ref.fasta.fai -Sb Rep1_from_apo_c.sam | express --output-dir rep1_express_f1_from_apo_c modwithindel_ref.fasta - &
samtools view -t modwithindel_ref.fasta.fai -Sb Rep1_from_ir64_c.sam | express --output-dir rep1_express_f1_from_ir64_c modwithindel_ref.fasta - &

samtools view -t modwithindel_ref.fasta.fai -Sb Rep1_from_apo_s.sam | express --output-dir rep1_express_f1_from_apo_s modwithindel_ref.fasta - &
samtools view -t modwithindel_ref.fasta.fai -Sb Rep1_from_ir64_s.sam | express --output-dir rep1_express_f1_from_ir64_s modwithindel_ref.fasta - &

samtools view -t modwithindel_ref.fasta.fai -Sb Rep2_from_apo_c.sam | express --output-dir rep2_express_f1_from_apo_c modwithindel_ref.fasta - &
samtools view -t modwithindel_ref.fasta.fai -Sb Rep2_from_ir64_c.sam | express --output-dir rep2_express_f1_from_ir64_c modwithindel_ref.fasta - &

samtools view -t modwithindel_ref.fasta.fai -Sb Rep2_from_apo_s.sam | express --output-dir rep2_express_f1_from_apo_s modwithindel_ref.fasta - &
samtools view -t modwithindel_ref.fasta.fai -Sb Rep2_from_ir64_s.sam | express --output-dir rep2_express_f1_from_ir64_s modwithindel_ref.fasta - &

