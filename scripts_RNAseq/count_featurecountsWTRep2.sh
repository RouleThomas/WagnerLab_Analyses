#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=30G



nthreads=12



input_list=("WT_Rep2")

for x in "${input_list[@]}"; do
    # run featurecounts

featureCounts -p -C -O -a ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chr.gtf -o counts/${x}.txt mapped_STAR_trim/${x}Aligned.sortedByCoord.out.bam
done