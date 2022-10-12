# Re-analyses of the Tan et al 2022 RNAseq using:
[Paper](https://academic.oup.com/plcell/article-abstract/34/8/2969/6580212?redirectedFrom=fulltext&login=false#367433007)

- trimming: Trimmomatic 
- mapping: STAR
- counting: featurecounts
- DE: Deseq2


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


## PCA on bam files ##
To check how well the data cluster. Notably, what is this WT_Rep3, good to use?\
1. Need create a bed file containing all exons
Already perform with hisat2 mapping *../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/exons.tsv*
2. Use multiBamSummary (deeptools) to generate a matrix of read counts
Install deeptools in CondaGS environment and launch command.\
```bash
module load Anaconda/2019.10
conda activate CondaGS
pip install deeptools
multiBamSummary BED-file --BED ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/exons.tsv --bamfiles mapped_STAR/SDG711RNAi_Rep1Aligned.sortedByCoord.out.bam mapped_STAR/SDG711RNAi_Rep2Aligned.sortedByCoord.out.bam mapped_STAR/WT_Rep1Aligned.sortedByCoord.out.bam mapped_STAR/WT_Rep3Aligned.sortedByCoord.out.bam mapped_STAR_trim/WT_Rep2Aligned.sortedByCoord.out.bam -o SDG_WT_matrix.npz
```
Launch as ```sbatch scripts/matrix_bam.sh```; Submitted batch job 200425=DONE, 
3. Use PlotPCA to make PCA plot
Example:\
```bash
plotPCA -in [multiBamSummary_prefix].npz \
    --transpose \
    --ntop 0 \
    --labels [samp_1] [samp_2] [...] [samp_n] \
    -o [plotPCA_prefix].png \
    -T "PCA of read counts in exons"
```
- ```transpose``` Transpose matrix in the format of rows to be samples and the columns to be features
- ```--ntop 0``` All rows are used in PCA
- [samp] list should be respective to multiBamSummary `bamfiles` parameter
```bash
plotPCA -in SDG_WT_matrix.npz \
    --transpose \
    --ntop 0 \
    --labels mapped_STAR/SDG711RNAi_Rep1Aligned.sortedByCoord.out.bam mapped_STAR/SDG711RNAi_Rep2Aligned.sortedByCoord.out.bam mapped_STAR/WT_Rep1Aligned.sortedByCoord.out.bam mapped_STAR/WT_Rep3Aligned.sortedByCoord.out.bam mapped_STAR_trim/WT_Rep2Aligned.sortedByCoord.out.bam \
    -o plotPCA_SDG_WT.png \
    -T "PCA of read counts in gene"
```
Launched as ```sbatch scripts/plot_PCA.sh```; Submitted batch job 200580=DONE. Replicate two is weird, maybe explain why it has a different name since the begining...! 
![plot](https://github.com/RouleThomas/WagnerLab_Analyses/blob/main/data_RNAseq/plotPCA_SDG_WT.png)

**Conclusion:** Use WT Rep1 and Rep3 and RNAi Rep1 and Rep2 for Deseq2

## Generation of the coverage (wig) files ##
--> Generate coverage file (bigwig)from bam file
Aligned.sortedByCoord.out.bam files has not been indexed, script *mapping_STAR_raw.sh* has been corrected but indexation launch as ```sbatch scripts/index_bam.sh```; Submitted batch job 200418=DONE
```bash
bamCoverage --bam mapped_STAR/${x}Aligned.sortedByCoord.out.bam --outFileName mapped_STAR/${x}.Aligned.sortedByCoord.out.bw --outFileFormat bigwig --normalizeUsing BPM --binSize 10  
```
- --binSize 10 will give a 10bp resolution, file may be too big, if that is the case, lets increase it to 50bp
- normalization in BPM=TPM
Always better to normalize per TPM. Launch as ```sbatch scripts/BamToBigwig.sh```; Submitted batch job 200421=DONE, files looks good.


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
- Relaunch hisat2 mapping on raw files to count with it see wether STAR is the isse ```scripts/mapping_hisat2_raw.sh```; Submitted batch job 198231=DONE.
- Relaunch STAR mapping on raw files using the good gtf to see if using different GTF for mapping and counting may be the isse; 1st re-index the genome ```sbatch scripts/STAR_indexation.sh```; Submitted batch job 198233=FAIL, fasta and gtf chr name not the same (fasta genome name is chr01 vs gtf is 1...) 
**troubleshooting solution:** Make fasta and GTF file same header... Lets put FASTA header as 1,2,etc... So here is command that 1st remove *chr0* and then remove *chr* only remaining for the chr>10: ```sed 's/chr0//' ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_genome.fasta | sed 's/chr//' > ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_genome_numeric.fasta```; ```sbatch scripts/STAR_indexation.sh```, Submitted batch job 198234=DONE
Here command to add *chr* in front of the 1,2,3,etc of the gtf (but useless in my case) ```awk '{ if($1 !~ /^#/){print "chr"$0} else{print $0} }' ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chr.gtf > ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chrlabel.gtf```\
Then re-mapping ```sbatch scripts/mapping_STAR_raw.sh```; Submitted batch job 198235=DONE, same fail (quality vs sequence lenght) for WT_Rep2. 
**troubleshooting solution:** Issue was that gtf annotation file use was different for counting and genome indexation, always use the same gtf file!!! Now that mapping is good (good gtf) compare counting between STAR and hisat2, test on *WT_Rep1*:
- ```featureCounts -p -C -O -P -B -d 30 -D 1000 -M --fraction```: multimapped reads counted STAR_raw: 20873031 77.8% 
- ```featureCounts -p -C -O -P -B -d 30 -D 1000```: multimapped reads not counted STAR_raw: 720043972 74.7%
- ```featureCounts -p -C -O```: multimapped reads not counted, non-paired reads counted, STAR_raw: 23685433 88.2%
- ```featureCounts -p -C -O```: multimapped reads not counted, non-paired reads counted, hisat2_raw: 3458904 12.0%
**Let's not count multimapped reads but count the non-paired reads using STAR.**,Launch as ```sbatch scripts/count_featurecounts.sh```, For WT_Rep2 trim_crop15bp has been used; Submitted batch job 200414=DONE. All good >85% Successfully assigned alignments.\
Except WT_Rep2 failed because come from mapping from another GTF.\
Need repeat Trimming (let's not crop, will gain more information) ```sbatch scripts/trimmomatic_noadapter_WTRep2.sh```; Submitted batch job 200415=DONE\
Need repeat Mapping ```sbatch scripts/mapping_STAR_trim_WTRep2.sh```; Submitted batch job 200416=DONE\
Re-counting ```sbatch scripts/count_featurecountsWTRep2.sh```; Submitted batch job 200419=DONE; 20819532 (88.0%).\
**All files are ready for the DEGs calculation**


## DEGs with Deseq2 ##
See [help](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). Input is unormalized fragment (as paired-end) counts. As I used featurecounts input is *Count matrix*.\
Let's set up an R environment for Deseq2. Let's use the *featurecounts* conda environment for Deseq2.\
1. Install the Deseq2 library using R script
```bash
Rscript scripts/install_R_Deseq2.R
```
Look like it failed, here is last rows of output:
```R
The downloaded source packages are in
        ‘/tmp/RtmpRNox3Y/downloaded_packages’
Installation paths not writeable, unable to update packages
  path: /cm/shared/apps/R/4.2.1/lib64/R/library
  packages:
    cluster, foreign, MASS, Matrix, nlme, nnet, survival
Warning messages:
1: In install.packages(...) :
  installation of package ‘RcppArmadillo’ had non-zero exit status
2: In install.packages(...) :
  installation of package ‘mixsqp’ had non-zero exit status
3: In install.packages(...) :
  installation of package ‘ashr’ had non-zero exit status
[1] "PACKAGE VERSIONS:"
[1] "argparse"
[1] ‘2.1.6’
[1] "tidyverse"
Error in packageVersion("tidyverse") :
  there is no package called ‘tidyverse’
Execution halted
```
Once in R, ```library("DESeq2")``` result in nothing....

**troubleshooting:** Lets not use the script to install, better go into R and install "manually", follow [bioconductor](https://bioconductor.org/packages/release/bioc/html/DESeq2.html):
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("DESeq2")
```
**troubleshooting:** Multiple fail, because of depedencies with other packages. Maybe also because I am on node04 in which I encountered issue to run CHIPQC in R in the past!\
Try install as follow a, within node01:
```R
BiocManager::install("DESeq2", dependencies = TRUE)
```
**troubleshooting:** ```Error: C++11 standard requested but CXX11 is not defined```\
modify the Makeconf file from R and add the following; 
- I do not know where is my R folder dependent on my conda environment (conda folder here : ```../../roule/.conda/envs/featurecounts/``` but no R inside...). 
- Try using a previous version of R manually installed where I see the Makeconf file. 
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("DESeq2")
```
Looks like it work!\
**troubleshooting:** R parameter version was not good, using an R version install locally within the cluster, it work. **so use ```/home/roule/R/R-4.2.0/bin/R``` to run DESeq2**
```R
library("DESeq2")
```
**troubleshooting solution:** weird character appear while deleting ```^H``` or using arrows... It is extremely annoying , so create a new conda environment using R4.2.0:
```bash
conda create --name DESeq2 r-base=4.2.0
conda activate DESeq2
```
Then force re-installation of DESeq2 ```BiocManager::install("DESeq2", force=TRUE)``` in R, ERROR: ``` there is no package called ‘Matrix’``` so I install it ```BiocManager::install("Matrix")```. Same with ```‘codetools’, ‘survival’``` and WORKS!!!\
**So to use DESeq2, use *DESeq2* conda environment**

2. Construct the DESeqDataSet using a count matrix (featurecounts output)
Tidy in R featurecount outputs, import file and keep only column gene ID and counts for each sample.\
```R
library("DESeq2")
library("tidyverse")
library("apeglm") # BiocManager::install("apeglm")
getwd()
SDG711RNAi_Rep1 = read_delim("counts/SDG711RNAi_Rep1.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
```
Create function to transform tibble into matrix:
```R
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
# Transform tibble into matrix
SDG711RNAi_Rep1_matrix = make_matrix(select(SDG711RNAi_Rep1, -Geneid), pull(SDG711RNAi_Rep1, Geneid))
```
Import all count sample (WT Rep 1 and 3), and combine into 1 file (row= gene and column= condition/replicate); then transform into a single matrix:
```R
# import and keep gene ID and counts
WT_Rep1 = read_delim("counts/WT_Rep1.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>% select(Geneid, `mapped_STAR/WT_Rep1Aligned.sortedByCoord.out.bam`) %>% rename("WT_Rep1"=`mapped_STAR/WT_Rep1Aligned.sortedByCoord.out.bam`)
WT_Rep3 = read_delim("counts/WT_Rep3.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>% select(Geneid, `mapped_STAR/WT_Rep3Aligned.sortedByCoord.out.bam`) %>% rename("WT_Rep3"=`mapped_STAR/WT_Rep3Aligned.sortedByCoord.out.bam`)
SDG711RNAi_Rep2 = read_delim("counts/SDG711RNAi_Rep2.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>% select(Geneid, `mapped_STAR/SDG711RNAi_Rep2Aligned.sortedByCoord.out.bam`) %>% rename("SDG711RNAi_Rep2"=`mapped_STAR/SDG711RNAi_Rep2Aligned.sortedByCoord.out.bam`)
SDG711RNAi_Rep1 = read_delim("counts/SDG711RNAi_Rep1.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>% select(Geneid, `mapped_STAR/SDG711RNAi_Rep1Aligned.sortedByCoord.out.bam`) %>% rename("SDG711RNAi_Rep1"=`mapped_STAR/SDG711RNAi_Rep1Aligned.sortedByCoord.out.bam`)

# merge all sample into one datafile
counts_all = WT_Rep1 %>% left_join(WT_Rep3) %>% left_join(SDG711RNAi_Rep1) %>% left_join(SDG711RNAi_Rep2)

# transform merge tibble into matrix
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 
```
Create data file ```coldata``` describing each samples; **columns of the count matrix and the rows of the column data in same order**:\
```R
# Generate sample description df
coldata_raw <- data.frame (sample  = c("WT_Rep1", "WT_Rep3", "SDG711RNAi_Rep1", "SDG711RNAi_Rep2"),
                  genotype = c("WT", "WT", "SDG77RNAi", "SDG77RNAi") )

# transform df into matrix
coldata = make_matrix(select(coldata_raw, -sample), pull(coldata_raw, sample))

# Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct
```
Construct the DESeqDataSet:
```R
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design = ~ genotype)
```
3. Perform the DEGs

```R
# Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

# Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

# Differential expression analyses
dds <- DESeq(dds)
res <- results(dds) # No need to use 'contrast' as only 1 comparison (eg. WT vs mutant)
res 

# Export result as 'raw_table'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="DESeq2/raw_table.txt")
raw_table <- read_csv("C:/Users/roule/Box/Ongoing projects/Cluster_GPC/Tan et al 2022 RNAseq/data_RNAseq/DESeq2/raw_table.txt") #To import

# Shrunken LFC for better estimation (notably for low counts and high disperson)
resLFC <- lfcShrink(dds, coef="genotype_SDG77RNAi_vs_WT", type="apeglm")

# Export result as 'shrunken_table'
write.csv(resLFC %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="DESeq2/shrunken_table.txt")

# Some summary (numb of DEGs)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
```
Output statistics
```
out of 32137 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 1576, 4.9%
LFC < 0 (down)     : 812, 2.5%
outliers [1]       : 0, 0%
low counts [2]     : 1870, 5.8%
(mean count < 4)
2388 DEGs
```
Can also do Independent hypothesis weighting, another method to obtain adjusted pvalue, more stringent tho... (do it when time, need troubleshoot ```BiocManager::install("IHW")``` installation!)


4. Data vizualization
*plotMA* LFC variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1.\
Better to do it on shrunken log2 fold changes (as it remove noise associated with LFC from low count genes).
```R
plotMA(resLFC, ylim=c(-2,2))
```
**troubleshooting:** Looks like X11 is disabled using ```capabilities()```\
```
Error in .External2(C_X11, d$display, d$width, d$height, d$pointsize,  :
  unable to start device X11cairo
In addition: Warning message:
In (function (display = "", width, height, pointsize, gamma, bg,  :
  unable to open connection to X11 display ''
```

Lets **back up all files** to come back later after troubleshooting X11:
```R
write.table(res,file="DESeq2/test") # save as matrix
write.table(resLFC,file="DESeq2/resLFC") # save as matrix
write.table(counts_all_matrix,file="DESeq2/counts_all_matrix") # save as matrix
res = read.table("DESeq2/res",header=TRUE,row.names=1) # load matrix
resLFC = read.table("DESeq2/resLFC",header=TRUE,row.names=1) # load matrix
counts_all_matrix = read.table("DESeq2/counts_all_matrix",header=TRUE,row.names=1) # load matrix
```


*PCA plot*

There is many other plot, explore [doc](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input) to know more.







