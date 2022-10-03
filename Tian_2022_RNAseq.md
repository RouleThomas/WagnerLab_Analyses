# Re-analyses of the Tan et al 2022 RNAseq using:
[Paper](https://academic.oup.com/plcell/article-abstract/34/8/2969/6580212?redirectedFrom=fulltext&login=false#367433007)

- trimming: Trimmomatic 
- mapping: STAR
- counting: XXX
- DE: XXX


## 1. import raw fastq file ##
--> Setup prerequisets and Greenscreen conda environment
```bash
srun --mem=20g --pty bash -l
cd Tian_2022TPC_ChIP
module load sratoolkit/2.11.2
module load Anaconda/2019.10
conda activate CondaGS
```

 
--> Download files from the SRA (paired end mode):\
*Files information can be found in RNAseq analyses Tan et al 2022.xlsx*
```bash
fasterq-dump SRR15663633 -S
fasterq-dump SRR15663632 -S
fasterq-dump SRR15663619 -S
fasterq-dump SRR15663618 -S
fasterq-dump SRR15663623 -S
```
ongoing

## 2. Start analyses ##
--> Compress/tidy/rename file with .sh script (script edited with nano)
```bash
sbatch scripts/organize_raw_fasta.sh
sbatch scripts/organize_raw_fasta1.sh
```
Submitted batch job 198123=DONE\
Submitted batch job 198124=DONE

--> Check presence of adaptors with FASTQC\
*Kmer Content and Sequence Duplication levels  always fail for RNA-seq*\
```bash
sbatch scripts/fastqc_raw.sh
```
Submitted batch job 198129=DONE




```bash
sbatch scripts/fastqc_raw.sh
```
Submitted batch job 198129=

###Test both Trimming and non-Trimming of the reads, then follow with the method showing best mapping###
--> Index the genome with STAR
```bash
module load STAR/2.7.10a
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ../GreenScreen/rice/GreenscreenProject/meta/genome_STAR_RNAseq --genomeFastaFiles ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_genome.fasta --sjdbGTFfile ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/transcripts_exon.gff 
```
Comand launch into ```sbatch scripts/STAR_indexation.sh```\
Submitted batch job 198132=DONE, got the following warning ```!!!!! WARNING: --genomeSAindexNbases 14 is too large for the genome size=373245519, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 13```
So re-run with the correction indicated
Submitted batch job 198133=CHUI AL

## --> Trimming reads and mapping ##




## --> Non-trimming of the reads and mapping ##
Go for mapping directly with the unmapped raw reads:\
Command use at AVIESAN SCHOOL:
```bash
STAR --genomeDir /shared/bank/eba2018/star/Arabidopsis_thaliana.TAIR10.1 --runThreadN 4 --readFilesCommand zcat --outFileNamePrefix output_STAR/${SAMPLE}_ --readFilesIn ${INPUT} --outSAMtype BAM SortedByCoordinate --alignIntronMax 1000
```
Command used by Sammy:
```bash
STAR --readFilesIn /home/wagner-lab/sklasfeld/Projects/TFL1/RNAseq/V2/trimmed_fastq/ft_24hrFRP_R1.trimmed.fastq --outFileNamePrefix ft_24hrFRP_R1 --runThreadN 12 --runMode alignReads --genomeDir /home/wagner-lab/sklasfeld/Araport11/STAR_genome_dir --sjdbOverhang 84 --outSAMprimaryFlag AllBestScore --outSJfilterCountTotalMin 10 5 5 5 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --alignIntronMin 60 --alignIntronMax 6000 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 2 --outWigType bedgraph
```
We need to adapt as we are working with Rice genome\













