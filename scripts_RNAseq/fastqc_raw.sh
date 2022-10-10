#!/bin/bash

# create output directory
mkdir -p data/FASTQC/raw



x=("SDG711RNAi_Rep1_1" "SDG711RNAi_Rep2_1" "WT_Rep1_1" "WT_Rep2_1" "WT_Rep3_1" "SDG711RNAi_Rep1_2" "SDG711RNAi_Rep2_2" "WT_Rep1_2" "WT_Rep2_2" "WT_Rep3_2")
        

for x in "${x[@]}"; do
    ../GreenScreen/Software/fastqc_v0.11.9/FastQC/fastqc -o data/FASTQC/raw fastq/raw/${x}.fastq.gz
done