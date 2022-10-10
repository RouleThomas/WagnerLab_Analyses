#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=20G
#SBATCH --nodelist=node01



input_list=("EMF2_Rep1" "EMF2_Rep2" "H3K27me3_Rep1" "H3K27me3_Rep2")

for x in "${input_list[@]}"; do
                                                 
bamCoverage --bam mapped/chip/downsample/${x}.dupmark.sorted.bam --outFileName mapped/chip/downsample/${x}.dupmark.sorted.bw --outFileFormat bigwig --binSize 1 --numberOfProcessors 4 --extendReads --scaleFactor 0.5                                          
                                                 
done                                  