#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=40G


# hard-code directories

fastqc_trim_dir="data/FASTQC/trimmed"



# run fastqc oldschool way


    ../GreenScreen/Software/fastqc_v0.11.9/FastQC/fastqc -o ${fastqc_trim_dir} fastq/trimmed/SDG711RNAi_Rep1_1.trimmed_paired.fastq.gz 
    ../GreenScreen/Software/fastqc_v0.11.9/FastQC/fastqc -o ${fastqc_trim_dir} fastq/trimmed/SDG711RNAi_Rep1_2.trimmed_paired.fastq.gz 
    ../GreenScreen/Software/fastqc_v0.11.9/FastQC/fastqc -o ${fastqc_trim_dir} fastq/trimmed/SDG711RNAi_Rep2_1.trimmed_paired.fastq.gz 
    ../GreenScreen/Software/fastqc_v0.11.9/FastQC/fastqc -o ${fastqc_trim_dir} fastq/trimmed/SDG711RNAi_Rep2_2.trimmed_paired.fastq.gz 

    ../GreenScreen/Software/fastqc_v0.11.9/FastQC/fastqc -o ${fastqc_trim_dir} fastq/trimmed/WT_Rep1_1.trimmed_paired.fastq.gz 
    ../GreenScreen/Software/fastqc_v0.11.9/FastQC/fastqc -o ${fastqc_trim_dir} fastq/trimmed/WT_Rep1_2.trimmed_paired.fastq.gz 
    