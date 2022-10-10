#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=40G



nthreads=12
module load samtools/1.15


# create output directory
mkdir -p mapped_STAR

input_list=("SDG711RNAi_Rep1" "SDG711RNAi_Rep2" "WT_Rep1" "WT_Rep3")

for x in "${input_list[@]}"; do
    # run STAR
samtools index mapped_STAR/${x}Aligned.sortedByCoord.out.bam 
done
