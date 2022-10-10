#!/bin/bash

# create output directory
mkdir -p data/FASTQC/raw


x=("EMF2_Rep1_1" "EMF2_Rep1_2" "EMF2_Rep2_1" "EMF2_Rep2_2" "H3K27me3_Rep1_1" 
   "H3K27me3_Rep1_2" "H3K27me3_Rep2_1" "H3K27me3_Rep2_2" "IgG_1" "IgG_2"
   "input_1" "input_2")
        

for x in "${x[@]}"; do
    ../GreenScreen/Software/fastqc_v0.11.9/FastQC/fastqc -o data/FASTQC/raw fastq/raw/${x}.fastq.gz
done
