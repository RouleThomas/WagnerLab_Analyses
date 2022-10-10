#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=30G



nthreads=4
module load STAR/2.7.10a


# create output directory
mkdir -p mapped_STAR_trimmed_crop15bp

input_list=("SDG711RNAi_Rep1" "SDG711RNAi_Rep2" "WT_Rep1" "WT_Rep2" "WT_Rep3")

for x in "${input_list[@]}"; do
    # run STAR
STAR --genomeDir ../GreenScreen/rice/GreenscreenProject/meta/genome_STAR_RNAseq --runThreadN 4 --readFilesCommand zcat --outFileNamePrefix mapped_STAR_trimmed_crop15bp/${x} --readFilesIn fastq/trimmed_crop15bp/${x}_1.trimmed_paired.fastq.gz fastq/trimmed_crop15bp/${x}_2.trimmed_paired.fastq.gz --outSAMtype BAM SortedByCoordinate
done

