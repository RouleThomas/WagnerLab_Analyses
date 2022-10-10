#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=30G



nthreads=4
module load samtools/1.15


# create output directory
mkdir -p mapped_hisat2

input_list=("SDG711RNAi_Rep1" "SDG711RNAi_Rep2" "WT_Rep1" "WT_Rep2" "WT_Rep3")

for x in "${input_list[@]}"; do
    # run hisat2 mapping
../GreenScreen/Software/hisat2/hisat2 -x ../GreenScreen/rice/GreenscreenProject/meta/genome_hisat2_RNAseq -1 fastq/raw/${x}_1.fastq.gz -2 fastq/raw/${x}_2.fastq.gz -S mapped_hisat2/${x}.sam
    # sort the reads
samtools sort -o mapped_hisat2/${x}.bam \
        mapped_hisat2/${x}.sam
    # index the bam file
samtools index mapped_hisat2/${x}.bam 
done


