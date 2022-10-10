#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=30G



input_list=("SDG711RNAi_Rep1" "SDG711RNAi_Rep2" "WT_Rep1" "WT_Rep3")

for x in "${input_list[@]}"; do
                                                 
bamCoverage --bam mapped_STAR/${x}Aligned.sortedByCoord.out.bam --outFileName mapped_STAR/${x}.Aligned.sortedByCoord.out.bw --outFileFormat bigwig --normalizeUsing BPM --binSize 10                             
done                                  

bamCoverage --bam mapped_STAR_trim/WT_Rep2Aligned.sortedByCoord.out.bam --outFileName mapped_STAR_trim/WT_Rep2.Aligned.sortedByCoord.out.bw --outFileFormat bigwig --normalizeUsing BPM --binSize 10   