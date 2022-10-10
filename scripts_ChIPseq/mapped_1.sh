#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=20G
#SBATCH --nodelist=node03


nthreads=3
PICARD_PATH="../GreenScreen/Software/picard/build/libs"
GENOME_PATH="../GreenScreen/rice/GreenscreenProject/meta/genome/bowtie2_genome_dir/IRGSP"
module load samtools/1.15
module load bowtie2/2.4.5



# create output directory
mkdir -p mapped/chip

input_list=("EMF2_Rep1" "EMF2_Rep2")


for x in "${input_list[@]}"; do
    # run bowtie2
    bowtie2  --phred33 -q \
	-x ${GENOME_PATH} \
        -S mapped/chip/${x}.sam \
        -1 fastq/trimmed/${x}_1.trimmed_paired.fastq.gz \
        -2 fastq/trimmed/${x}_2.trimmed_paired.fastq.gz
    # sort the reads
    samtools sort -o mapped/chip/${x}.bam \
        mapped/chip/${x}.sam
    # index the bam file
    samtools index mapped/chip/${x}.bam
    # remove reads without MAPQ>=30
    samtools view -@ ${nthreads} -F 772 -q 30 \
        -b mapped/chip/${x}.bam \
	chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 | \
        samtools sort -o mapped/chip/${x}.filter.bam
    # index filtered reads
    samtools index mapped/chip/${x}.filter.bam
    # mark duplicates with picard
    java -jar ${PICARD_PATH}/picard.jar MarkDuplicates \
        -I mapped/chip/${x}.filter.bam \
        -O mapped/chip/${x}.dupmark.bam \
        -M mapped/chip/${x}.dup.qc \
        -VALIDATION_STRINGENCY LENIENT \
        -REMOVE_DUPLICATES false \
	-ASSUME_SORTED true
    # sort reads after marking the duplicates
    samtools sort -o mapped/chip/${x}.dupmark.sorted.bam \
        mapped/chip/${x}.dupmark.bam
    # index the sorted reads
    samtools index mapped/chip/${x}.dupmark.sorted.bam
done
