#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=40G


nthreads=12

module load STAR/2.7.10a

STAR --runThreadN 12 --runMode genomeGenerate --genomeSAindexNbases 13 \
--genomeDir ../GreenScreen/rice/GreenscreenProject/meta/genome_STAR_RNAseq \
--genomeFastaFiles ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_genome_numeric.fasta \
--sjdbGTFfile ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chr.gtf 
