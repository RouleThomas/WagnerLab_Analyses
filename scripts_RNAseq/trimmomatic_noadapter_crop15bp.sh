#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=40G


# hard-code directories
fastq_raw_dir="fastq/raw"
fastq_trim_dir="fastq/trimmed_crop15bp"
fastqc_trim_dir="data/FASTQC/trimmed_crop15bp"

# create output directories
mkdir -p ${fastq_trim_dir}
mkdir -p ${fastqc_trim_dir}

trimmomatic_install_dir="../GreenScreen/Software"

# run Trimmomatic without adapter removal and crop 15bp extremity

no_adapter_input_list=("SDG711RNAi_Rep1" "SDG711RNAi_Rep2" "WT_Rep1" "WT_Rep2" "WT_Rep3") 

for x in "${no_adapter_input_list[@]}"; do
    raw_fastq1="${fastq_raw_dir}/${x}_1.fastq.gz"
    raw_fastq2="${fastq_raw_dir}/${x}_2.fastq.gz"
    trim_fastq1="${fastq_trim_dir}/${x}_1.trimmed_paired.fastq.gz"
    trim_fastq2="${fastq_trim_dir}/${x}_1.trimmed_unpaired.fastq.gz"
    trim_fastq3="${fastq_trim_dir}/${x}_2.trimmed_paired.fastq.gz"
    trim_fastq4="${fastq_trim_dir}/${x}_2.trimmed_unpaired.fastq.gz"
    java -jar ${trimmomatic_install_dir}/Trimmomatic-0.39/Trimmomatic-0.39/trimmomatic-0.39.jar \
        PE -threads 3 -phred33 ${raw_fastq1} ${raw_fastq2} \
        ${trim_fastq1} ${trim_fastq2} ${trim_fastq3} ${trim_fastq4}\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33 HEADCROP:15

    ${trimmomatic_install_dir}/fastqc_v0.11.9/FastQC/fastqc -o ${fastqc_trim_dir} ${trim_fastq1}
    ${trimmomatic_install_dir}/fastqc_v0.11.9/FastQC/fastqc -o ${fastqc_trim_dir} ${trim_fastq2}
    ${trimmomatic_install_dir}/fastqc_v0.11.9/FastQC/fastqc -o ${fastqc_trim_dir} ${trim_fastq3}
    ${trimmomatic_install_dir}/fastqc_v0.11.9/FastQC/fastqc -o ${fastqc_trim_dir} ${trim_fastq4}
done