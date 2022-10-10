#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=40G



nthreads=12


multiBamSummary BED-file --BED ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/exons.tsv --bamfiles mapped_STAR/SDG711RNAi_Rep1Aligned.sortedByCoord.out.bam mapped_STAR/SDG711RNAi_Rep2Aligned.sortedByCoord.out.bam mapped_STAR/WT_Rep1Aligned.sortedByCoord.out.bam mapped_STAR/WT_Rep3Aligned.sortedByCoord.out.bam mapped_STAR_trim/WT_Rep2Aligned.sortedByCoord.out.bam -o SDG_WT_matrix.npz