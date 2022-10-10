#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=30G



nthreads=4
module load STAR/2.7.10a


# create output directory
mkdir -p mapped_STAR_trim

input_list=("WT_Rep2")

for x in "${input_list[@]}"; do
    # run STAR
STAR --genomeDir ../GreenScreen/rice/GreenscreenProject/meta/genome_STAR_RNAseq --runThreadN 12 --readFilesCommand zcat --outFileNamePrefix mapped_STAR_trim/${x} --readFilesIn fastq/trimmed/${x}_1.trimmed_paired.fastq.gz fastq/trimmed/${x}_2.trimmed_paired.fastq.gz --outSAMtype BAM SortedByCoordinate
done