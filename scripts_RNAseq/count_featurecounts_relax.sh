#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=30G



nthreads=12



input_list=("SDG711RNAi_Rep1" "SDG711RNAi_Rep2" "WT_Rep1" "WT_Rep2" "WT_Rep3")

for x in "${input_list[@]}"; do
    # run featurecounts
featureCounts -p -C -O -a ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chr.gtf -o counts_relax/${x}.txt mapped_STAR/${x}Aligned.sortedByCoord.out.bam
done
