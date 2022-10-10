#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=40G


nthreads=12


../GreenScreen/Software/hisat2/hisat2-build ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_genome.fasta ../GreenScreen/rice/GreenscreenProject/meta/genome_hisat2_RNAseq --ss ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/splicesites.tsv --exon ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/exons.tsv
 