#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=40G



nthreads=12
module load STAR/2.7.10a
module load samtools/1.15

# create output directory
mkdir -p mapped_STAR

input_list=("SDG711RNAi_Rep1" "SDG711RNAi_Rep2" "WT_Rep1" "WT_Rep2" "WT_Rep3")

for x in "${input_list[@]}"; do
    # run STAR
STAR --genomeDir ../GreenScreen/rice/GreenscreenProject/meta/genome_STAR_RNAseq --runThreadN 12 --readFilesCommand zcat --outFileNamePrefix mapped_STAR/${x} --readFilesIn fastq/raw/${x}_1.fastq.gz fastq/raw/${x}_2.fastq.gz --outSAMtype BAM SortedByCoordinate
samtools index mapped_STAR/${x}Aligned.sortedByCoord.out.bam 
done
