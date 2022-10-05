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
*Kmer Content and Sequence Duplication levels always fail for RNA-seq*\
```bash
sbatch scripts/fastqc_raw.sh
```
Submitted batch job 198129=DONE\
FAIL: Per base sequence content/Per sequence GC content/Sequence Duplication Levels
-  Per base sequence content: The random hexamer primers, which are used to generate the cDNA library from your RNA transcripts were shown to not bind completly random. This non-random binding leads to this bias in "per base sequence content" from base 1-15.
-  Sequence Duplication Levels: can be because of highle expressed genes
-  Per sequence GC content: Can be adaptors, trimming may remove that, or contamination (that will not mapped to genome in that case)
Quality looks good for direct mapping without mapping -> no presence of adaptor so let's use standard trimming parameter


###Test both Trimming and non-Trimming of the reads, then follow with the method showing best mapping###
--> Index the genome with STAR
```bash
module load STAR/2.7.10a
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ../GreenScreen/rice/GreenscreenProject/meta/genome_STAR_RNAseq --genomeFastaFiles ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_genome.fasta --sjdbGTFfile ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/transcripts_exon.gff 
```
Comand launch into ```sbatch scripts/STAR_indexation.sh```\
Submitted batch job 198132=DONE, got the following warning ```!!!!! WARNING: --genomeSAindexNbases 14 is too large for the genome size=373245519, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 13```
So script re-run with the correction indicated\
Submitted batch job 198133=DONE (no warning)

## Trimming reads and mapping ##
-->Trimming\
Lets try to use the standard parameter as for ChIP\
```bash
java -jar ${trimmomatic_install_dir}/Trimmomatic-0.39/Trimmomatic-0.39/trimmomatic-0.39.jar \
        PE -threads 3 -phred33 ${raw_fastq1} ${raw_fastq2} \
        ${trim_fastq1} ${trim_fastq2} ${trim_fastq3} ${trim_fastq4}\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33
```
Launch as ```sbatch scripts/trimmomatic_noadapter.sh``` and ```sbatch scripts/trimmomatic_noadapter1.sh```
Submitted batch job 198153=DONE\
Submitted batch job 198154=DONE\
Trimming performed succesfully, fastqc failed, let's re-run a script for fastqc ```sbatch scripts/trimmomatic_fastqc.sh``` only while correcting the initial script ```scripts/trimmomatic_noadapter.sh``` for future use
Submitted batch job 198162=DONE\
Submitted batch job 198163=DONE\
**Trimming does not change anything (same fastqc FAIL)**. Maybe I did not trim enough, lets change parameter so that it trimm out 15nt (as the 10 first from *Per base sequence content* looks bad): added ```HEADCROP:15``` parameter:
```bash
java -jar ${trimmomatic_install_dir}/Trimmomatic-0.39/Trimmomatic-0.39/trimmomatic-0.39.jar \
        PE -threads 3 -phred33 ${raw_fastq1} ${raw_fastq2} \
        ${trim_fastq1} ${trim_fastq2} ${trim_fastq3} ${trim_fastq4}\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33 HEADCROP:15
```
Launched as ```sbatch scripts/trimmomatic_noadapter_crop15bp.sh```; Submitted batch job 198165\
trimming OK, some fastqc fail (failed: SDG711RNAi_Rep1_1.trimmed_paired.fastq.gz, SDG711RNAi_Rep1_1.trimmed_unpaired.fastq.gz, Failed to process file SDG711RNAi_Rep1_2.trimmed_paired.fastq.gz, SDG711RNAi_Rep1_2.trimmed_unpaired.fastq.gz, SDG711RNAi_Rep2_2.trimmed_paired.fastq.gz, SDG711RNAi_Rep2_2.trimmed_unpaired.fastq.gz, WT_Rep1_1.trimmed_unpaired.fastq.gz, WT_Rep1_2.trimmed_unpaired.fastq.gz, WT_Rep2_2.trimmed_paired.fastq.gz, WT_Rep2_2.trimmed_unpaired.fastq.gz, WT_Rep3_1.trimmed_paired.fastq.gz, WT_Rep3_1.trimmed_unpaired.fastq.gz, WT_Rep3_2.trimmed_paired.fastq.gz, WT_Rep3_2.trimmed_unpaired.fastq.gz ```Failed to process file SDG711RNAi_Rep1_1.trimmed_paired.fastq.gz
uk.ac.babraham.FastQC.Sequence.SequenceFormatException: Ran out of data in the middle of a fastq entry.  Your file is probably truncated
	at uk.ac.babraham.FastQC.Sequence.FastQFile.readNext(FastQFile.java:179)
	at uk.ac.babraham.FastQC.Sequence.FastQFile.next(FastQFile.java:125)
	at uk.ac.babraham.FastQC.Analysis.AnalysisRunner.run(AnalysisRunner.java:77)
	at java.lang.Thread.run(Thread.java:748)```
 **troubleshooting:** looks like I got fastqc for all so should be ok for interpretation\
 **Trimming has solved the previous *Per base sequence content* FAILED** But there is less reads so not sure that is the best approaches (raw read may be better, lets see!)
 
In parrallel, test trimming using PE adapter file:
```bash
java -jar ${trimmomatic_install_dir}/Trimmomatic-0.39/Trimmomatic-0.39/trimmomatic-0.39.jar \
        PE -threads 3 -phred33 ${raw_fastq1} ${raw_fastq2} \
        ${trim_fastq1} ${trim_fastq2} ${trim_fastq3} ${trim_fastq4}\
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33
```

Launched as ```sbatch scripts/trimmomatic_PEadapter.sh```; Submitted batch job 198166=DONE\
Looks similar as without using the adapters.\
**Let's mapp the 15bp crop reads** with STAR and hisat2:
- STAR default parameter (trim reads) (```sbatch scripts/mapping_STAR_trimmed.sh```; Submitted batch job 198186=CHUI AL





## Non-trimming of the reads and mapping ##
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
Let s use the mapping parameter from the [Science paper](https://www.science.org/doi/full/10.1126/science.aax8862) from Julia Bailey-Serres.\
They used Bowtie2/hisat2 allowing 2 nt mismatch (Tan et al paper use the same tool). --> Lets use Hisat2 (as bowtie2 do not support spliced alignemnt)


- STAR default parameter (raw reads) (```sbatch scripts/mapping_STAR_raw.sh```; Submitted batch job 198164=DONE
```STAR --genomeDir ../GreenScreen/rice/GreenscreenProject/meta/genome_STAR_RNAseq --runThreadN 4 --readFilesCommand zcat --outFileNamePrefix mapped_STAR/${x} --readFilesIn fastq/raw/${x}_1.fastq.gz fastq/raw/${x}_2.fastq.gz --outSAMtype BAM SortedByCoordinate```
**success except one file**: ```EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length SOLUTION: fix your fastq file``` for WT_rep2, Lets try to investigate the fastq file!

CHUI AL

- hisat2 default parameter (raw reads)
Installation following [hisat2 github](https://github.com/DaehwanKimLab/hisat2).\
Genome indexation launch under ```scripts/hisat2_indexation.sh```; Submitted batch job 198167=CHUI AL
















