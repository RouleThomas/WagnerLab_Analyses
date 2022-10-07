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
- STAR default parameter (trim reads) (```sbatch scripts/mapping_STAR_trimmed.sh```; Submitted batch job 198186=DONE, here no error message regarding fastq from WT_Rep2





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
**success except one file**: ```EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length SOLUTION: fix your fastq file``` for WT_rep2, Lets try to investigate the fastq file! Could have been due to the loop according to [forum](https://github.com/alexdobin/STAR/issues/1055)! lets try to remap that file only, same command; Submitted batch job 198195=SAME FAIL. Lets look at the fastq file where error detected: ```zgrep -A4 "@SRR15663632.24281996" fastq/raw/WT_Rep2_1.fastq.gz```, looks like the quality string lenght is equal to sequence lenght...
**troubleshooting**:I tried re-run the command within the cluster not as a slurm job and obtain same error... Lets try to unzip fastq and check for format error at this specific line and remove using ```nano``` and ```CTRL+W``` [ID line](https://github.com/alexdobin/STAR/issues/726). FAIL never try to open a fastq with nano :D.\
**troubleshooting**:Lets try to do the mapping with the uncompressed fastq file:
Submitted batch job 198204=FAIL again
**troubleshooting**: Check if the fastq is corrupted using ```cat fastq/raw/WT_Rep2_1.fastq | paste - - - - | awk -F '\t' '(length($2)!=length($4))'``` replace cat by ```gunzip -c *.fq.gz``` if working with compressed fastq: If output not empty, file is corrupted (lenght data quality not the same as lenght sequence). No output... The file is not corrupted!!!
**troubleshooting solution**: let's use the trim file for this file


- hisat2 default parameter (raw reads)
Installation following [hisat2 github](https://github.com/DaehwanKimLab/hisat2).\
Genome indexation launch under ```scripts/hisat2_indexation.sh```; Submitted batch job 198167=DONE\
```bash
../GreenScreen/Software/hisat2/hisat2 -x ../GreenScreen/rice/GreenscreenProject/meta/genome_hisat2_RNAseq -1 fastq/raw/${x}_1.fastq.gz -2 fastq/raw/${x}_2.fastq.gz -S mapped_hisat2/${x}.sam
```
**Let's redo it and provide the gff while building the index, that should improve the mapping...**
```
../GreenScreen/Software/hisat2/hisat2-build ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_genome.fasta ../GreenScreen/rice/GreenscreenProject/meta/genome_hisat2_RNAseq --ss --exon ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/transcripts_exon.gff
```
Submitted batch job 198198=FAIL ```Error: could not open --exon``` 
**troubleshooting**: maybe need to convert GFF to GTF, lets do it with [gffread](https://github.com/gpertea/gffread)\
```bash
../GreenScreen/Software/gffread/gffread ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/transcripts_exon.gff -T -o ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/transcripts_exon.gtf
```
Re run indexation with the gtf ```sbatch scripts/hisat2_indexation.sh```; Submitted batch job 198199=FAIL\
**troubleshooting**: issue is I need a splice site file and an exon file for --ss and --exon, and not GTF. [I can make these files using python script from hisat2](https://rnabio.org/module-01-inputs/0001/04/01/Indexing/):\
```bash
../GreenScreen/Software/hisat2/hisat2_extract_splice_sites.py ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/transcripts_exon.gtf > ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/splicesites.tsv
../GreenScreen/Software/hisat2/hisat2_extract_exons.py ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/transcripts_exon.gtf > ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/exons.tsv
```
FAIL, result in empty file...\
**troubleshooting**: Probably my gtf is not format as hisat2-extract wants... Lets look at the python script more in detail
Lets try with Arabidopsis gff, see if same fail. 
```
../GreenScreen/Software/gffread/gffread ../GreenScreen/tutorial/GreenscreenProject/meta/ArabidopsisGenome/Araport11_GFF3_genes_transposons.201606.gff -T -o ../GreenScreen/tutorial/GreenscreenProject/meta/ArabidopsisGenome/Araport11_GFF3_genes_transposons.201606.gtf
../GreenScreen/tutorial/GreenscreenProject/meta/ArabidopsisGenome/exons.tsv
```
That work! So the gtf file is indeed misformated in rice. Looks like in each row of rice gtf gene-id is not indicated. The Gff looks like shit also. Lets download another one: *Oryza_sativa.IRGSP-1.0.54.chr.gff3.gz* from [here](https://plants.ensembl.org/Oryza_sativa/Info/Index).\
```bash
../GreenScreen/Software/gffread/gffread ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chr.gff3 -T -o ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chr.gtf
../GreenScreen/Software/hisat2/hisat2_extract_splice_sites.py ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chr.gtf > ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/splicesites.tsv
../GreenScreen/Software/hisat2/hisat2_extract_exons.py ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chr.gtf > ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/exons.tsv
```
Work, proceed with the indexation
```bash
sbatch scripts/hisat2_indexation.sh
```
Submitted batch job 198201=DONE
**troubleshooting conclusion:** The gff file was not well formated! Lets use the Oryza_sativa.IRGSP-1.0.54.chr.gtf from now on for hisat2 (the gff is weird as it shows the chromosome as a gene...). Otherwise "transcripts_exons.gff" used for STAR mapping is exactly the same.\

--> Mapping with index genome specifying exons and splice sites\
Launch for raw and trimmed reads: ```sbatch scripts/mapping_hisat2_raw.sh``` *Submitted batch job 198222* and ```sbatch scripts/mapping_hisat2_crop15bp.sh``` *Submitted batch job 198214*, also convertion into bam file and indexed them.\

**Mapping conclusion comparison tool and parameter:**
- STAR (check the *Log.final.out file*): raw reads show more uniquely mapped reads than trim reads!
- hisat2 (no output so take value from the slurm log): raw reads are better. But much fewer uniquely mapped reads as compare to STAR raw.
Lets perform the counting using STAR raw mapping. All other mapped files deleted to save space...


## Generation of the coverage (wig) files ##



## Counting ##
Let's check if the data are stranded or non-stranded looking at a bam file on IGV. Color alignment per read strand on IGV showed specific orientation so **data is unstranded**.\

Sammy used ht-seqcount, I used featureCounts, and Julia used summarizeOverlaps function from the GenomicRanges package (weird!!!!) and Tan XXX.\
From documentation, featureCount is more adapted for paired-end data. Indeed if a reads map on two different gene, it will be ambiguous in htseq count, even if two of the paired reads map on one of the two gene!\
Lets use a conda environment for featurecount:\
```bash
module load Anaconda/2019.10
conda create -c bioconda -n featurecounts subread
conda activate featurecounts
```
Parameter to use:
- paired-end mode -p (add -C to not count paired reads with each pair on different chromosome!)
- -O to assign reads to metafeatures, so genes
- -P -B -d 30 -D 1000 (set min and max paired reads to 30 to 1000bp)
- stranded? No
- -M --fraction to count multimapped reads fractionally? Not recommended for RNAseq so I do not do it. In addition most of my reads are uniquely mapped
```bash
featureCounts -p -C -O -P -B -d 30 -D 1000 -a ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chr.gtf -o counts/${x}.txt mapped_STAR/${x}Aligned.sortedByCoord.out.bam
```
Launched as ```sbatch scripts/count_featurecounts.sh```; Submitted batch job 198225=cancel\
Also test without ```-P -B``` parameters, to count non both end mapped reads, thus more reads count:\
``` bash
featureCounts -p -C -O -a ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chr.gtf -o counts_relax/${x}.txt mapped_STAR/${x}Aligned.sortedByCoord.out.bam
```
Launched as ```sbatch scripts/count_featurecounts_relax.sh```; Submitted batch job 198227=cancel\
Show very low (<25%) of *Successfully assigned alignments* for both jobs... Looks like majority of reads do not align due to no features... 
**troubleshooting:** Lets try to allow multimapped reads adding ```-M --fraction``` parameter; Submitted batch job 198230=cancel, same low percent...
**troubleshooting:** 
- Relaunch hisat2 mapping on raw files to count with it see wether STAR is the isse ```scripts/mapping_hisat2_raw.sh```; Submitted batch job 198231=XXX
- Relaunch STAR mapping on raw files using the good gtf to see if using different GTF for mapping and counting may be the isse; 1st re-index the genome ```sbatch scripts/STAR_indexation.sh```; Submitted batch job 198233=FAIL, fasta and gtf chr name not the same (fasta genome name is chr01 vs gtf is 1...) 
**troubleshooting solution:** Make fasta and GTF file same header... Lets put FASTA header as 1,2,etc... So here is command that 1st remove *chr0* and then remove *chr* only remaining for the chr>10: ```sed 's/chr0//' ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_genome.fasta | sed 's/chr//' > ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_genome_numeric.fasta```; ```sbatch scripts/STAR_indexation.sh```, Submitted batch job 198234=DONE
Here command to add *chr* in front of the 1,2,3,etc of the gtf (but useless in my case) ```awk '{ if($1 !~ /^#/){print "chr"$0} else{print $0} }' ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chr.gtf > ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chrlabel.gtf```\
Then re-mapping ```sbatch scripts/mapping_STAR_raw.sh```; Submitted batch job 198235=XXX

Count the new hisat2 and STAR mapping, test with and without multimapping count and if shit again, maybe that is ribosomal RNA remaining???



























