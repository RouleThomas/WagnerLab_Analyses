#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=30G



nthreads=4
module load STAR/2.7.10a


# create output directory
mkdir -p mapped_STAR

input_list=("WT_Rep2")

for x in "${input_list[@]}"; do
    # run STAR
STAR --genomeDir ../GreenScreen/rice/GreenscreenProject/meta/genome_STAR_RNAseq --runThreadN 4 --outFileNamePrefix mapped_STAR/${x} --readFilesIn fastq/raw/${x}_1.fastq fastq/raw/${x}_2.fastq --outSAMtype BAM SortedByCoordinate
done