# Re-analyses of the Tian et al 2022 ChIP using the greenscreen pipeline
- Tan et al 2022 paper [here](https://academic.oup.com/plcell/article/34/8/2969/6580212#367433007)
- Klasfeld et al 2022 paper [here](https://www.biorxiv.org/content/10.1101/2022.02.27.482177v1)\
Tian et al used IgG as control for EMF2b ChIP, let's try using input instead

- trimming: Trimmomatic 
- mapping: bowtie2
- peak calling: MACS2
- assin peak to gene (custom python script)

## 1. import raw fastq file ##
--> Setup prerequisets and Greenscreen conda environment
```bash
srun --mem=20g --pty bash -l
cd Tian_2022TPC_ChIP
module load sratoolkit/2.11.2
module load Anaconda/2019.10
conda activate CondaGS
```
Installed sound effect to be notify on command runnin more than 5seconds by following [this](https://github.com/c0rp-aubakirov/notify_after_command_executed/)
```bash
 git clone https://github.com/c0rp-aubakirov/notify_after_command_executed.git
 cd notify_after_command_executed/
 echo "source $(pwd)/postexec_notify" >> ~/.bashrc
 ```
 Restart terminal and tested if notification appear after 5sec ```sleep 6``` : \
 @FAIL:  Does not work when I restart, mention "postexec_notify: No such file or directory"
 
--> Download files from the SRA (paired end mode):\
*Files information can be found in ChIP analyses Tan et al 2022.xlsx*
```bash
fasterq-dump SRR18596327 -S
fasterq-dump SRR18596325 -S
fasterq-dump SRR15663624 -S
fasterq-dump SRR15663617 -S
fasterq-dump SRR15663616 -S
fasterq-dump SRR15663629 -S
```
## 2. Start analyses (Peak calling broad H3K27me3 and Narrow EMF2) ##
*20220915*\
--> Compress/tidy/rename file with .sh script (script edited with ```nano```)
```bash
sbatch scripts/organize_raw_fasta.sh
```
Submitted batch job 197628=DONE\
*20220916*\
--> FASTQC to check adaptors type
```bash
sbatch scripts/fastqc_raw.sh
```
Submitted batch job 197643=DONE\
--> Trimming (some have adapters; some not; see *Tan et al 2022.xlsx*)
```bash
sbatch scripts/trimmomatic.sh
sbatch scripts/trimmomatic_adaptor_removal.sh
```
Submitted batch job 197661=DONE   \
Submitted batch job 197662=DONE   \

*20220919* \
Over-representated sequence still present in H3K27me3 IP \
Blast it on NCBI(**AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA**TCATTCCAGACAACTCCT): bold region blast to SARS-Cov2.... \
*I think ok to keep these, they will just not align to the genome*\

--> Mapping
```bash
sbatch scripts/mapped_1.sh
```
Submitted batch job 197710 (FIE)=DONE\
Submitted batch job 197712 (H3K27me3)=DONE\
Submitted batch job 197714 (IGG)=DONE\
Submitted batch job 197715 (input)=DONE\
Start/End 1.20pm\10pm End \ 

--> Downsample replicate (for H3K27me3 and EMF2 2 Replicates)
```bash
sbatch scripts/downsampleBam.sh
```
Submitted batch job 197750\
Start/End 11.48am/=FAIL\
**output:**\
*[SEVERE][Biostar145820]Problem writing temporary file file:///tmp/sortingcollection.3596660749531177790.tmp.  Try setting TMP_DIR to a file system with lots of space.
htsjdk.samtools.util.RuntimeIOException: Problem writing temporary file file:///tmp/sortingcollection.3596660749531177790.tmp.  Try setting TMP_DIR to a file system with lots of space.*\
**troubleshoots:**\
- Try create a tmp folder in working directory to store temporary files and edit java comand (follow [this](https://www.biostars.org/p/42613/)) 
```bash
mkdir tmp
java -Djava.io.tmpdir=`pwd`/tmp TMP_DIR=`pwd`/tmp -jar COMAND
```
Submitted batch job 197762=FAIL\
**output:**\
*Error: Could not find or load main class TMP_DIR=.home.roule.Tian_2022TPC_ChIP.tmp*\
**troubleshoots:**\
- Try direct to tmp folder directly (follow [this](https://www.biostars.org/p/42613/)) 
```bash
java -Djava.io.tmpdir=tmp TMP_DIR=tmp -jar COMAND
```
Submitted batch job 197763=FAIL\
**output:**\
Error: Could not find or load main class TMP_DIR=tmp
**troubleshoots:**\
- Try to put the tmp-related argument at the end of the command
Submitted batch job 197764=WORK BUT NEW FAIL\
**output:**\
```
[SEVERE][Biostar145820]There was an error in the command line arguments. Please check your parameters : Expected one or zero argument but got 3 : [mapped/chip/EMF2_Rep2.dupmark.sorted.bam, -Djava.io.tmpdir=tmp, TMP_DIR=tmp]
com.github.lindenb.jvarkit.lang.JvarkitException$CommandLineError: There was an error in the command line arguments. Please check your parameters : Expected one or zero argument but got 3
```
**troubleshoots:**\
- Try to re-order correctly the different arguments and add ```-Xmx2g``` to allow 2G max of temporary files
```bash
 java -Xmx2g -Djava.io.tmpdir=tmp -jar ${JVARKIT_PATH}/dist/biostar145820.jar \
     --seed ${seed} -n ${min_val} \
     -o ${out_dir}/${samp}_Rep${rep}.dupmark.bam \
     ${in_dir}/${samp}_Rep${rep}.dupmark.sorted.bam \
     TMP_DIR=tmp
```
Submitted batch job 197767=FAIL\
**output:**\
```
[SEVERE][Biostar145820]There was an error in the command line arguments. Please check your parameters : Expected one or zero argument but got 2 : [mapped/chip/EMF2_Rep2.dupmark.sorted.bam, TMP_DIR=tmp]
```
**troubleshoots:**\
- Remove the ```TMP_DIR=tmp```
Submitted batch job 197769=DONE\
Last=2hrs for 2 sample with 2 replicates\
**Troubleshootings Conclusion:** The java syntax was not good; For me ```-Xmx2g -Djava.io.tmpdir=tmp``` (store temporary file to tmp folder and allow 2G max) needs to be in between ```java``` (call java) and ```-jar``` (JAR file is used: biostar145820.jar); Then other arguments ```-n``` minimum number of read in one of the replicate; ```--seed 42``` (because that is the number of the universe, lol, could be any number... That is to generate randomness and reproducibility as a subset of the reads are taken)

--> Call peak IP to input control individually\
```bash
module load Anaconda/2019.10
conda activate CondaGS
sbatch scripts/macs2_callpeaks.sh
```
Submitted batch job 197779=FAIL\
**output:**\
```
INFO  @ Tue, 20 Sep 2022 17:43:30: #3 Call peaks for each chromosome... 
ValueError: cannot resize this array: it does not own its data
```
**troubleshoots:**\
Seems OK to ignore the error. Nevertheless, seems to be cause by not enough space in temporary shared folder, let's try adding one with the option ```--tempdir PATH```
```bash
module load Anaconda/2019.10
conda activate CondaGS
sbatch scripts/macs2_callpeaks_test.sh
```
Submitted batch job 197801=DONE but other FAIL\
No more error message but files are the same so that could have been ignore.\
**output:**\
*awk: `10' argument to `-v' not in 'var=value' form*
Need to add ```q=${q}``` after the -v argument, as follow: ```awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.narrowPeak > ${macs2_out}/noMask_qval${q}/${x}_peaks.narrowPeak``` >>> Script *macs2_callpeaks.sh* has been corrected

--> Select Broad peaks for histone and IgG\
Needs to modify the macs2 command (just added ```--broad``` parameter:
```bash
macs2_callpeaks_broad.sh
```
Submitted batch job 197820=FAIL at greenscreen filtering\
Submitted batch job 197831 to complete with greenscreen filtering=DONE\


--> Call peak (narrow and broad) on pooled bam
```bash
sbatch scripts/macs2_callpeaks_BamPool.sh
sbatch scripts/macs2_callpeaks_BamPool_broad.sh
```
Submitted batch job 197853 (narrow) =DONE 
Submitted batch job 197854 (broad) =DONE


--> Generate coverage file (bigwig)from bam file\
Try [bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) since we are using Paired-end sequencing data\
```bash
bamCoverage \
--bam <input>.bam \
--outFileName <output>.bw \
--outFileFormat bigwig \
--binSize 1 \
--numberOfProcessors 4 \
--extendReads
--scaleFactor 0.5
```
```--binSize 1``` for good resolution; ```--scaleFactor 0.5``` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1.; ```--extendReads``` Reads extented taking into account mean fragment size of all mated reads.
```bash
sbatch scripts/BamToBigwig.sh
sbatch scripts/BamToBigwig1.sh
```
Submitted batch job 197810 (EMF2, H3K27me3)=DONE\
Submitted batch job 197811 (input, IgG)=DONE\

--> Merge the two bigwig replicates into 1 bigwig
wiggletool merge parameter, install it from git:
```bash
git clone https://github.com/dpryan79/libBigWig.git
cd libBigWig
make install
git clone --recurse-submodules https://github.com/samtools/htslib.git
cd htslib 
make install
```
FAIL at make for each command, authorization stuff\
Try install through Conda instead:\
Go to new conda environment just in case\
```bash
conda activate CondaUmap
conda install -c bioconda wiggletools
```
DONE, bigwig mean > bedGraph > bigwig
```bash
wiggletools mean mapped/chip/downsample/H3K27me3_Rep1.dupmark.sorted.bw mapped/chip/downsample/H3K27me3_Rep2.dupmark.sorted.bw > mapped/chip/downsample/H3K27me3.bedGraph
```
FAIL, format is not in bedgraph : write_bg
```bash
wiggletools write_bg mapped/chip/downsample/H3K27me3.wig mean mapped/chip/downsample/H3K27me3_Rep1.dupmark.sorted.bw mapped/chip/downsample/H3K27me3_Rep2.dupmark.sorted.bw
/home/roule/GreenScreen/Software/bedGraphToBigWig mapped/chip/downsample/H3K27me3.wig ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_chr_count.txt mapped/chip/downsample/H3K27me3.bw
```
DONE\
**Troubleshootings Conclusion:** Command to use to pool 2 bigwig:
```bash
conda activate CondaUmap
wiggletools write_bg mapped/chip/downsample/H3K27me3.wig mean mapped/chip/downsample/H3K27me3_Rep1.dupmark.sorted.bw mapped/chip/downsample/H3K27me3_Rep2.dupmark.sorted.bw
/home/roule/GreenScreen/Software/bedGraphToBigWig mapped/chip/downsample/H3K27me3.wig ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_chr_count.txt mapped/chip/downsample/H3K27me3.bw
```
**Alternative:** I may also make a wig from the pooled bam file

--> Venn Diagram mine and Tan analyses
Import the bed2venn.py python script to generate venn diagram from two bed file from [here](https://github.com/YichaoOU/HemTools/)
```bash
python3 scripts/bed2venn.py -b1 data/peaks_for_comparison/H3K27me3_pool_peaks.broadPeak -b2 data/peaks_for_comparison/H3K27me3_peaks_Tan.bed -l1 Greenscreen -l2 Tan et al
```
**output:**\
*File "scripts/bed2venn.py", line 88
    print "A: Number of lines in %s: %s"%(args.b1,F1)
                                       ^
SyntaxError: invalid syntax*
**troubleshoots:**\
Change to python2 using Umap conda environment:
```bash
python scripts/bed2venn.py -b1 data/peaks_for_comparison/H3K27me3_pool_peaks.broadPeak -b2 data/peaks_for_comparison/H3K27me3_peaks_Tan.bed -l1 Greenscreen -l2 Tan
```
**output:**\
```
ImportError: No module named matplotlib
ImportError: No module named seaborn
ImportError: No module named matplotlib_venn
```
**troubleshoots:**\
```bash
pip install matplotlib
pip install seaborn
pip install matplotlib_venn
```
Many errors:
```
ERROR: Received illegal bin number 4294967295 from getBin call.
ERROR: Unable to add record to tree.
Traceback (most recent call last):
  File "scripts/bed2venn.py", line 105, in <module>
    main()
  File "scripts/bed2venn.py", line 87, in main
    out = wccount(outFile)
  File "scripts/bed2venn.py", line 25, in wccount
    df = pd.read_csv(filename,sep="\t",header=None)
  File "/home/roule/.conda/envs/CondaUmap/lib/python2.7/site-packages/pandas/io/parsers.py", line 702, in parser_f
    return _read(filepath_or_buffer, kwds)
  File "/home/roule/.conda/envs/CondaUmap/lib/python2.7/site-packages/pandas/io/parsers.py", line 429, in _read
    parser = TextFileReader(filepath_or_buffer, **kwds)
  File "/home/roule/.conda/envs/CondaUmap/lib/python2.7/site-packages/pandas/io/parsers.py", line 895, in __init__
    self._make_engine(self.engine)
  File "/home/roule/.conda/envs/CondaUmap/lib/python2.7/site-packages/pandas/io/parsers.py", line 1122, in _make_engine
    self._engine = CParserWrapper(self.f, **self.options)
  File "/home/roule/.conda/envs/CondaUmap/lib/python2.7/site-packages/pandas/io/parsers.py", line 1853, in __init__
    self._reader = parsers.TextReader(src, **kwds)
  File "pandas/_libs/parsers.pyx", line 545, in pandas._libs.parsers.TextReader.__cinit__
pandas.errors.EmptyDataError: No columns to parse from file*
```
**troubleshoots:**\
- Try correct/clean the input beds and keep the first three column of both bed files
```bash
cat data/peaks_for_comparison/H3K27me3_peaks_Tan.bed | awk '{ if ($2!=$3) print $0 }' > data/peaks_for_comparison/H3K27me3_peaks_Tan_corr.bed
cat data/peaks_for_comparison/H3K27me3_pool_peaks.broadPeak | awk '{ if ($2!=$3) print $0 }' > data/peaks_for_comparison/H3K27me3_pool_corr.bed
```
- also tried:
```bash
awk '($2<$3){print $0}' data/peaks_for_comparison/H3K27me3_peaks_Tan.bed > data/peaks_for_comparison/H3K27me3_peaks_Tan_corr.bed
awk '($2<$3){print $0}' data/peaks_for_comparison/H3K27me3_peaks_Tan > data/peaks_for_comparison/H3K27me3_peaks_Tan_corr.bed
```
**ouputs:**\
New ERROR:
```
ModuleCmd_Load.c(213):ERROR:105: Unable to locate a modulefile for 'bedtools'
```
It fail at running bedtools at line 82, perform slight modification (use bedops instead of bedtools) and that run but very few overlap but they should be a lot! I think because bedops --intersect not doing exactly same stuff as bedtools intersect -u
**troubleshoots:**\
Let's try [this](https://github.com/asntech/intervene) instead:\
Install through Conda to have all dependency; create a new conda environment for Venn:\
```bash
conda create --name venn
conda activate venn
conda install -c bioconda intervene
intervene --help
```
**ouputs:**\
pkg_resources.ContextualVersionConflict: (numpy 1.14.2 (/home/roule/.conda/envs/venn/lib/python3.6/site-packages), Requirement.parse('numpy>=1.15'), {'seaborn'})
**troubleshoots:**\
```bash
pip3 install 'numpy==1.15.0' #install numpy1.15 and unistalled numpy1.14 
```
DONE\
Launch the Venn command:\
```bash
intervene venn -i data/peaks_for_comparison/H3K27me3_pool_corr.bed data/peaks_for_comparison/H3K27me3_peaks_Tan_corr.bed --output data/peaks_for_comparison/H3K27me3
intervene venn -i data/peaks_for_comparison/EMF2_pool_peaks.narrowPeak data/peaks_for_comparison/EMF2_peaks_Tan.bed --output data/peaks_for_comparison/EMF2
intervene venn -i data/peaks_for_comparison/H3K27me3_pool_corr.bed data/peaks_for_comparison/EMF2_pool_peaks.narrowPeak --output data/peaks_for_comparison/EMF2_H3K27_greenscreen
```
DONE, files are in pdf format
Tested other representation\
```bash
intervene upset -i data/peaks_for_comparison/H3K27me3_pool_corr.bed data/peaks_for_comparison/H3K27me3_peaks_Tan_corr.bed data/peaks_for_comparison/EMF2_pool_peaks.narrowPeak data/peaks_for_comparison/EMF2_peaks_Tan.bed --output data/peaks_for_comparison/upset_EMF2_H3K27_all
```
**troubleshoots:**\
Needs to install manually into R UpsetR package
```R
install.packages("UpSetR")
q()
```
Launch the previously generated R script to generate the plot:\
```bash
/home/roule/R/R-4.2.0/bin/Rscript data/peaks_for_comparison/upset_EMF2_H3K27_all/Intervene_upset.R
```
DONE, generated plot in pdf

## Conclusion 1 ##

- H3K27me3 peaks identified are quite comparable in between greenscreen and Tan pipeline (XX% overlap)
- Much more EMF2 peaks identified using greenscreen pipeline as compare to Tan analyses (XX% more)
--> Let's now remove the H3K27me3 and EMF2 peaks that overlap with IgG peaks see what we obtain


## 3. Pursue analyses (Peak calling broad IgG and filtered out overlapping H3K27me3 and EMF2 peaks) ##
Call peaks IgG to input\
```bash
sbatch scripts/macs2_callpeaks_broad_IGG.sh
sbatch scripts/macs2_callpeaks_broad1.sh
```
Submitted batch job 197875= DONE, FAIL at filtering p-value and removing greenscreen region (corrected in the script macs2_callpeaks_broad_IGG.sh, but launch aditional to correct)\
Submitted batch job 197876= DONE\

Only 28 peaks has been called! That is very few, thus likely not change downstream analyses...\
To make it change, let's decrease the qvalue for calling the peaks from 10 to 5\
```bash
sbatch scripts/macs2_callpeaks_broad_IGG.sh # q value edited from 10 to 5 
```
Submitted batch job 197926=DONE\
Similarly, only 53 peaks identified.\
--> let's only remove the H3K27me3 and EMF2b peaks overlapping with the 28 IgG peaks using bedtools:
```bash
bedtools intersect -a data/macs2_out/chipPeaks/broad_gsMask_qval10/H3K27me3_pool_peaks.broadPeak -b data/macs2_out/chipPeaks/broad_gsMask_qval10/IgG_peaks.broadPeak -v > data/peaks_for_comparison/H3K27me3_noIgG.broadPeak
```



## 4. Annotate peaks to genes ##
Let's use [ChIPseeker](https://academic.oup.com/bioinformatics/article/31/14/2382/255379) to assign peak to genes:\
Lets create a specific ChIPseeker conda environment, to make sure I do not mess up my DESeq2 environment...\
1. Create a conda environment with all R capabilities TRUE; For that let's copy our DESeq2 conda environment specificity where everythings fine 

```bash
conda activate DESeq2
conda list --explicit > spec-file-DESeq2.txt # That is all our spec from the current conda env (here DESeq2)
conda create --name ChIPseeker --file spec-file-DESeq2.txt # Here create environment based on DESeq2 env spec files
```
**`conda list --explicit > spec-file-DESeq2.txt` is great command to backup environment**


2. Install ChIPseeker package in R
```R
BiocManager::install("ChIPseeker")
```
**troubleshooting solution:** `Configuration failed because libudunits2.so was not found`; lets install it using `conda install udunits2`.
**troubleshooting solution:** `configure: error: sf is not compatible with GDAL versions below 2.0.1`; lets install it using `conda install r-sf`.
```R
library("ChIPseeker")
library("tidyverse")
library("GenomicRanges") # BiocManager::install("GenomicRanges")
```
WORK!!\
Import and tidy ChIP peaks (narrow for EMF2, broad for H3K27me3):
```R
# Import
peaks_EMF2 =  read.table('data/macs2_out/chipPeaks/gsMask_qval10/EMF2_pool_peaks.narrowPeak') %>% rename(Chr=V1, start=V2, end=V3, name=V4, score=V5, strand=V6, signal_value=V7, pvalue=V8, qvalue=V9, peak=V10) # Import and rename columns

# Tidy peaks
peaks.gr = makeGRangesFromDataFrame(peaks_EMF2,keep.extra.columns=TRUE)

# Quick stat (peak per chr)
table(as.character(seqnames(peaks.gr)))

chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12
 2017  1535  1643  1224  1163  1097  1106  1011   862   783   878   801
 
 # plot (coverage distribution over genome)
 pdf('data/ChIPseeker/EMF2_coverage.pdf')
 covplot(peaks.gr,weightCol='fold_enrichment')
 dev.off()
 ```
 
**troubleshoot:** Trouble to make the txdb format, tried from either GTF or grange file result in the same error\
Need convert GTF into txdb format [this](https://rdrr.io/bioc/GenomicFeatures/man/makeTxDbFromGRanges.html) R function.
```R
# prepare packages
library('tracklayer')
library('checkr') # install.packages("checkr")
source("scripts/makeTxDbFromGRanges.R") # load 'makeTxDbFromGRanges' function
source("scripts/makeTxDb.R") # # load 'makeTxDb' function

# import gtf and transform into grange file
gr <- import("../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chr.gtf")

# make txdb from grange
txdb <- makeTxDbFromGRanges(gr)

# make txdb from grange only transcript metafeature
gr <- import("../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chr.gtf") # import gtf as grange file
gr_transcript = gr[ gr$type == "transcript" ] # keep only transcript feature
txdb <- makeTxDbFromGRanges(gr_transcript)

#make txdb from ensembl
library('RMariaDB') # BiocManager::install("RMariaDB")
source("scripts/Ensembl-utils.R")
library(RCurl) # install.packages("RCurl")


# mke txdb from gtf directly 
txdb <- makeTxDbFromGFF("../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_representative/Oryza_sativa.IRGSP-1.0.54.chrlabel.gtf", format= "gtf")
```
FAIL: Same error message obtain: `Make the TxDb object ... Error: exclusive must be a flag`; cannot find anything on internet...
Check how look a `txdb` format using `BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")`
 
**troubleshoot solution:** Generate txdb from Biomart db work! see [here](https://support.bioconductor.org/p/100772/)
 ```R
 # Import the txdb
 txdb <- makeTxDbFromBiomart(biomart="plants_mart",
                            dataset="osativa_eg_gene",
                            host="plants.ensembl.org")
                            
 # define promoter regions as +/- 3kb/500bp
 promoter = getPromoters(TxDb=txdb, upstream=3000, downstream=500)

## get the reads around the promoter regions
tagMatrix = getTagMatrix(peaks.gr, windows=promoter)
 ```
 **troubleshoot solution:** chromosome name for `peaks.gr` (*chr01*) and `promoter` (*1*) are not the same, change promoter into 1. Change the `peaks.gr` file accordingly
 ```R
# import 
peaks_EMF2 =  read.table('data/macs2_out/chipPeaks/gsMask_qval10/EMF2_pool_peaks.narrowPeak') %>% rename(Chr=V1, start=V2, end=V3, name=V4, score=V5, strand=V6, signal_value=V7, pvalue=V8, qvalue=V9, peak=V10)
 
#dataframe for correct chromosome name
chr_label <- data.frame (Chr  = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"),
                  chr = c(1:12)
                  )
                  
# join and select
peaks_EMF2_chrValues = peaks_EMF2 %>% left_join(chr_label) %>% dplyr::select(-Chr) %>% dplyr::select(chr, everything()) # join both dataframe, remove the previous bad Chr label and put the new chr label as first column

peaks.gr = makeGRangesFromDataFrame(peaks_EMF2_chrValues,keep.extra.columns=TRUE)
 
# get the reads around the promoter regions
tagMatrix = getTagMatrix(peaks.gr, windows=promoter)

# plot this as density heatmap
pdf('data/ChIPseeker/EMF2_heatmap.pdf')
tagHeatmap(tagMatrix, xlim=c(-3000, 500), color="red")
dev.off()

# plot this as profile plot
pdf('data/ChIPseeker/EMF2_plotAvgProf.pdf')
plotAvgProf(tagMatrix, xlim=c(-3000, 500), conf=0.95,resample=500, facet="row") # confidence interval is estimated by bootstrap method (500 iterations)
dev.off()

# Annotate peak to genes
peakAnno  = annotatePeak(peaks.gr,tssRegion=c(-3000,500), TxDb=txdb)
peak.anno = as.data.frame(peakAnno)
write.table(peak.anno,file='data/ChIPseeker/EMF2_peaks.anno.csv',row.names=FALSE,quote=FALSE,sep='\t') # Save table

# Check distribution plot relative to features
pdf('data/ChIPseeker/EMF2_distribution.pdf')
plotAnnoPie(peakAnno)
dev.off()

# Check distribution to gene
pdf('data/ChIPseeker/EMF2_distribution_gene.pdf')
plotPeakProf2(peak = peaks.gr, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 100,
              TxDb = txdb, ignore_strand = F) # nbin = filter out peak smaller than 100
dev.off()
```

Now lets do the same for **H3K27me3 broad** peaks:
 ```R
# import 
peaks_H3K27me3 =  read.table('data/macs2_out/chipPeaks/broad_gsMask_qval10/H3K27me3_pool_peaks.broadPeak') %>% rename(Chr=V1, start=V2, end=V3, name=V4, score=V5, strand=V6, signal_value=V7, pvalue=V8, qvalue=V9) # no column 10 peak when call --broad

 
#dataframe for correct chromosome name
chr_label <- data.frame (Chr  = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"),
                  chr = c(1:12)
                  )
                  
# join and select
peaks_H3K27me3_chrValues = peaks_H3K27me3 %>% left_join(chr_label) %>% dplyr::select(-Chr) %>% dplyr::select(chr, everything()) # join both dataframe, remove the previous bad Chr label and put the new chr label as first column

peaks.gr = makeGRangesFromDataFrame(peaks_H3K27me3_chrValues,keep.extra.columns=TRUE)
 
# get the reads around the promoter regions
tagMatrix = getTagMatrix(peaks.gr, windows=promoter)

# plot this as profile plot
pdf('data/ChIPseeker/H3K27me3_plotAvgProf.pdf')
plotAvgProf(tagMatrix, xlim=c(-3000, 500), conf=0.95,resample=500, facet="row") # confidence interval is estimated by bootstrap method (500 iterations)
dev.off()

# Annotate peak to genes
peakAnno  = annotatePeak(peaks.gr,tssRegion=c(-3000,500), TxDb=txdb)
peak.anno = as.data.frame(peakAnno)
write.table(peak.anno,file='data/ChIPseeker/H3K27me3_peaks.anno.csv',row.names=FALSE,quote=FALSE,sep='\t') # Save table

# Check distribution plot relative to features
pdf('data/ChIPseeker/H3K27me3_distribution.pdf')
plotAnnoPie(peakAnno)
dev.off()

# Check distribution to gene
pdf('data/ChIPseeker/H3K27me3_distribution_gene.pdf')
plotPeakProf2(peak = peaks.gr, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 100,
              TxDb = txdb, ignore_strand = F)
dev.off()
```
Lets  try to do venn diagram to see overlap EMF2 and H3K27me3 related to genes annotation follow [this](http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html)
```R
# Provide path as list for H3k27me3 and EMF2 peaks
files <- list(EMF2 = "data/macs2_out/chipPeaks/gsMask_qval10/EMF2_pool_peaks.narrowPeak", H3K27me3 = "data/macs2_out/chipPeaks/broad_gsMask_qval10/H3K27me3_pool_peaks.broadPeak")

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 500), verbose=FALSE)
```
**troubleshoot solution:** It fail as files and txdb not have similar seqname under grange... Lets correct
```R
## For EMF2 ##
# import 
peaks_EMF2 =  read.table('data/macs2_out/chipPeaks/gsMask_qval10/EMF2_pool_peaks.narrowPeak') %>% rename(Chr=V1, start=V2, end=V3, name=V4, score=V5, strand=V6, signal_value=V7, pvalue=V8, qvalue=V9, peak=V10)
 
#dataframe for correct chromosome name
chr_label <- data.frame (Chr  = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"),
                  chr = c(1:12)
                  )
                  
# join and select
peaks_EMF2_chrValues = peaks_EMF2 %>% left_join(chr_label) %>% dplyr::select(-Chr) %>% dplyr::select(chr, everything()) # join both dataframe, remove the previous bad Chr label and put the new chr label as first column

# save 
write.table(peaks_EMF2_chrValues,file='data/ChIPseeker/peaks_EMF2_chrValues',row.names=FALSE,quote=FALSE,sep='\t') # Save table


## For H3K27me3 ##
# import 
peaks_H3K27me3 =  read.table('data/macs2_out/chipPeaks/broad_gsMask_qval10/H3K27me3_pool_peaks.broadPeak') %>% rename(Chr=V1, start=V2, end=V3, name=V4, score=V5, strand=V6, signal_value=V7, pvalue=V8, qvalue=V9) # no column 10 peak when call --broad

 
#dataframe for correct chromosome name
chr_label <- data.frame (Chr  = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"),
                  chr = c(1:12)
                  )
                  
# join and select
peaks_H3K27me3_chrValues = peaks_H3K27me3 %>% left_join(chr_label) %>% dplyr::select(-Chr) %>% dplyr::select(chr, everything()) # join both dataframe, remove the previous bad Chr label and put the new chr label as first column

# save 
write.table(peaks_H3K27me3_chrValues,file='data/ChIPseeker/peaks_H3K27me3_chrValues',row.names=FALSE,quote=FALSE,sep='\t') # Save table

## Provide path as list for H3k27me3 and EMF2 peaks ##
files <- list(EMF2 = "data/ChIPseeker/peaks_EMF2_chrValues", H3K27me3 = "data/ChIPseeker/peaks_H3K27me3_chrValues")

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 500), verbose=FALSE)
```
We can make some graph for comparison:
```R
# plotAnnoBar
pdf('data/ChIPseeker/plotAnnoBar_comparison.pdf')
plotAnnoBar(peakAnnoList)
dev.off()

# plotDistToTSS
pdf('data/ChIPseeker/plotDistToTSS_comparison.pdf')
plotDistToTSS(peakAnnoList)
dev.off()
```
Check overlap between peaks and annotated genes
```R
genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)

# Venn diagram 
pdf('data/ChIPseeker/venn_peaks_genes_comparison.pdf')
vennplot(genes)
dev.off()

# Extract overlapping peaks
as_tibble(do.call(cbind, genes)) 
overlap_genes = intersect(genes$EMF2, genes$H3K27me3) %>% as_tibble() %>% unique() # Table of EMF2 and H3K27me3 bound genes

## Extract EMF2 peaks corresponding to these genes
# Import tidy
peaks_EMF2_chrValues = peaks_EMF2 %>% left_join(chr_label) %>% dplyr::select(-Chr) %>% dplyr::select(chr, everything()) # join both dataframe, remove the previous bad Chr label and put the new chr label as first column
peaks.gr = makeGRangesFromDataFrame(peaks_EMF2_chrValues,keep.extra.columns=TRUE)
peakAnno  = annotatePeak(peaks.gr,tssRegion=c(-3000,500), TxDb=txdb)
peak.anno = as_tibble(peakAnno)

# isolate EMF2 peak overlap genes and reorder as bed file with peak coordinate
overlap_genes_EMF2 = overlap_genes %>% rename(geneId=value) %>% left_join(peak.anno) %>% dplyr::select(seqnames, start, end, geneId, name, score, signal_value, pvalue, qvalue, peak, annotation, geneChr, geneStart, geneEnd, geneLength, geneStrand, transcriptId, distanceToTSS) %>% rename(chr=seqnames)

# Export
write.table(overlap_genes_EMF2, file="data/ChIPseeker/overlap_genes_EMF2_complete.bed",row.names=FALSE,quote=FALSE,sep='\t')
write.table(overlap_genes_EMF2 %>% dplyr::select(chr,start,end), file="data/ChIPseeker/overlap_genes_EMF2.bed",row.names=FALSE,quote=FALSE,sep='\t')
```
Check signficance overlap between peaks and annotated genes (in both direction)
```R
enrichPeakOverlap(queryPeak = "data/ChIPseeker/peaks_H3K27me3_chrValues", targetPeak = "data/ChIPseeker/peaks_EMF2_chrValues", nShuffle=1000, TxDb = txdb , pAdjustMethod="BH", chainFile=NULL)
enrichPeakOverlap(queryPeak = "data/ChIPseeker/peaks_EMF2_chrValues", targetPeak = "data/ChIPseeker/peaks_H3K27me3_chrValues", nShuffle=500, TxDb = txdb , pAdjustMethod="BH", chainFile=NULL)
```
Output (500 bootstrap last ~1h): 
Direction (targetPeak=EMF2); significant: `                   qSample              tSample qLen  tLen N_OL      pvalue
1 peaks_H3K27me3_chrValues peaks_EMF2_chrValues 5830 14120 5487 0.002364066
     p.adjust
1 0.002364066
`
Direction (targetPeak=H3K27me3):
`
1 peaks_EMF2_chrValues peaks_H3K27me3_chrValues 14120 5830 5487 0.002105263
     p.adjust
1 0.002105263
`


**Conclusion ChIPseeker:**
- Annotation peak to gene is done.
- Binding profile has been check (EMF2 is center to TSS whereas H3K27me3 is within gene body)
- Overlap EMF2 and H3K27me3 (related to peak in gene) has been check: significant overlap


## Motif discovery with MEME ##
Check [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/12_functional_analysis.html) and [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/motif_analysis_prep.html) for help\
Let's do **MEME-CHIP** (specifically design for ChIP, combine  DREME (motif discovery) and Tomtom (check if motif ressemble known TF))\
First need a FASTA containing sequence of our peaks. Let's focus on the EMF2 peaks from the genes that are enriched in H3K27me3 and DEGs\
Let's isolate these peaks in R:
```R
# import files
peaks_overlap =  read_delim('data/ChIPseeker/overlap_genes_EMF2_complete.bed')
peaks_overlap_tidy =  peaks_overlap %>% dplyr::select(chr, start, end, name, score, signal_value, qvalue, peak, geneId, distanceToTSS, annotation)

DEGs = read_delim('../Tian_2022TPC_RNAseq/DESeq2/raw_table.txt') %>% filter(padj<0.05) %>% dplyr::select(gene, baseMean, log2FoldChange, padj)
DEGs_tidy = DEGs %>% separate(gene, c("deleteme", "geneId"), sep=":") # remove "gene:" on gene column

# Combine DEGs with peaks overlap
DEGs_peakoverlap = DEGs_tidy %>% inner_join(peaks_overlap_tidy) %>% dplyr::select(chr, start, end, name, score, signal_value, qvalue, peak, distanceToTSS, annotation, geneId, baseMean, log2FoldChange, padj) # Tidy the order as bed

# Export the clean bed table (plus a version with chromsome name matching the FASTA file)
write.table(DEGs_peakoverlap, file="data/ChIPseeker/DEGs_peakoverlap.bed",row.names=FALSE,quote=FALSE,sep='\t')

chr_label <- data.frame (Chr  = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"),
                  chr = c(1:12)   )

DEGs_peakoverlap_tidy = DEGs_peakoverlap %>% left_join(chr_label) %>% dplyr::select(Chr, start, end, name, score, signal_value, qvalue, peak, distanceToTSS, annotation, geneId, baseMean, log2FoldChange, padj) # Tidy the order as bed
                  
write.table(DEGs_peakoverlap_tidy, file="data/ChIPseeker/DEGs_peakoverlap_tidy.bed",row.names=FALSE,quote=FALSE,sep='\t')
```

Vizualize peak distribution of these specific overlapping one, to compare with what I obtain broadly (EMF2 peaks solely) : XXX
```R
library(ChIPseeker)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- makeTxDbFromBiomart(biomart="plants_mart",
                            dataset="osativa_eg_gene",
                            host="plants.ensembl.org")
                            
# define promoter regions as +/- 3kb/500bp
promoter = getPromoters(TxDb=txdb, upstream=3000, downstream=500)

DEGs_peakoverlap =  read_delim('data/ChIPseeker/DEGs_peakoverlap.bed')


peaks.gr = makeGRangesFromDataFrame(DEGs_peakoverlap,keep.extra.columns=TRUE)

# Annotate peak to genes
peakAnno  = annotatePeak(peaks.gr,tssRegion=c(-3000,500), TxDb=txdb)

# Check distribution plot relative to features
pdf('data/ChIPseeker/DEGs_peakoverlap_distribution.pdf')
plotAnnoPie(peakAnno)
dev.off()

# Check distribution to gene
pdf('data/ChIPseeker/DEGs_peakoverlap_distribution_gene.pdf')
plotPeakProf2(peak = peaks.gr, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 100,
              TxDb = txdb, ignore_strand = F)
dev.off()
```
Usfell R function/package; biomaRt; can be found [here](https://github.com/grimbough/biomaRt)\

Convert bed to FASTA to be use in MEME-CHIP 
```bash
bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF>

bedtools getfasta -fi ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_genome.fasta -bed data/ChIPseeker/DEGs_peakoverlap_tidy_noheader.bed -fo data/ChIPseeker/DEGs_peakoverlap_tidy.fasta # BED header as been removed with manual editing
```
**MEME-CHIP:** The lenght of our DNA sequence is too long (505 DNA sequences, between 280 and 22325 in length (average length 2521.7)), should be 500bp to be optimal; so let's see wether size differ ~ features:
```R
library(tidyverse)
DEGs_peakoverlap =  read_delim('data/ChIPseeker/DEGs_peakoverlap.bed')

DEGs_peakoverlap_size = DEGs_peakoverlap %>% mutate(length = end-start) #obtain size of the peaks

# plot
ggplot(DEGs_peakoverlap_size %>% filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR","3' UTR", "Distal Intergenic", "Downstream (<=300bp)")), aes(annotation, length)) + geom_boxplot() + geom_hline(yintercept=500, linetype="dashed", color = "red", size=.5)

# Isolate small (400-800bp) PRE, specific to promoter region (from 0 to 5000bp to TSS)
DEGs_peakoverlap_size %>% filter(distanceToTSS %in% (0:5000), length %in% (400:800)) # n=34
DEGs_peakoverlap_size %>% filter(distanceToTSS %in% (0:5000), length < 1000) # n=55
DEGs_peakoverlap_size %>% filter(distanceToTSS %in% (0:5000), length < 1500) # n=104

# Export all tables to be converted in FASTA
## Tidy
DEGs_peakoverlap_size_400800 = DEGs_peakoverlap_size %>% filter(distanceToTSS %in% (0:5000), length %in% (400:800)) %>% left_join(chr_label) %>% dplyr::select(Chr, start, end, name, score, signal_value, qvalue, peak, distanceToTSS, annotation, geneId, baseMean, log2FoldChange, padj) # Tidy the order as bed
DEGs_peakoverlap_size_Less1000 = DEGs_peakoverlap_size %>% filter(distanceToTSS %in% (0:5000), length < 1000) %>% left_join(chr_label) %>% dplyr::select(Chr, start, end, name, score, signal_value, qvalue, peak, distanceToTSS, annotation, geneId, baseMean, log2FoldChange, padj) # Tidy the order as bed
DEGs_peakoverlap_size_Less1500 = DEGs_peakoverlap_size %>% filter(distanceToTSS %in% (0:5000), length < 1500) %>% left_join(chr_label) %>% dplyr::select(Chr, start, end, name, score, signal_value, qvalue, peak, distanceToTSS, annotation, geneId, baseMean, log2FoldChange, padj) # Tidy the order as bed

## Save
write.table(DEGs_peakoverlap_size_400800, file="data/ChIPseeker/DEGs_peakoverlap_size_400800.bed",row.names=FALSE,quote=FALSE,sep='\t')
write.table(DEGs_peakoverlap_size_Less1000, file="data/ChIPseeker/DEGs_peakoverlap_size_Less1000.bed",row.names=FALSE,quote=FALSE,sep='\t')
write.table(DEGs_peakoverlap_size_Less1500, file="data/ChIPseeker/DEGs_peakoverlap_size_Less1500.bed",row.names=FALSE,quote=FALSE,sep='\t')
```
Convert bed to FASTA to be use in MEME-CHIP 
```bash
bedtools getfasta -fi ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_genome.fasta -bed data/ChIPseeker/DEGs_peakoverlap_size_400800_noheader.bed -fo data/ChIPseeker/DEGs_peakoverlap_size_400800.fasta # BED header as been removed with manual editing
bedtools getfasta -fi ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_genome.fasta -bed data/ChIPseeker/DEGs_peakoverlap_size_Less1000_noheader.bed -fo data/ChIPseeker/DEGs_peakoverlap_size_Less1000_noheader.fasta
bedtools getfasta -fi ../GreenScreen/rice/GreenscreenProject/meta/genome/IRGSP-1.0_genome.fasta -bed data/ChIPseeker/DEGs_peakoverlap_size_Less1500_noheader.bed -fo data/ChIPseeker/DEGs_peakoverlap_size_Less1500_noheader.fasta
```
**MEME-CHIP:** Best output is for the 400-800 filter, even though n smaller. Shows AP2 TF family enrichment notably and other C2C2...

## Fine-tune motif finding ##
Modify my peak lenght: find the summit of each peaks and extend to +/- 300bp = 600bp total for each peaks; optimal for motif findings.\
Peak summit coordinates are found in the `data/macs2_out/chipPeaks/EMF2_pool_summits.bed`, let's inner_join summits from peaks corresponding to greenscreen called peaks `data/macs2_out/chipPeaks/gsMask_qval10/EMF2_pool_peaks.narrowPeak`; and extend to +/- 300bp in R:

```bash
srun --x11 --nodelist=node03 --mem=20g --pty bash -l
conda activate ChIPseeker
R
```
```R
library(ChIPseeker)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Open summit and greenscreen peaks files
EMF2_peaks_summit =  read_delim('data/macs2_out/chipPeaks/EMF2_pool_summits.bed', col_names=FALSE) %>% rename(chr=X1, start=X2, end=X3, name=X4, value=X5)
EMF2_peaks_greenscreen = read_delim('data/macs2_out/chipPeaks/gsMask_qval10/EMF2_pool_peaks.narrowPeak', col_names = FALSE) %>% rename(chr=X1, start=X2, end=X3, name=X4, value=X5) %>% dplyr::select(-X6, -X7, -X8, -X9, -X10)

# isolate the greenscreen EMF2 peaks from the summit files using peak name
EMF2_peaks_summit_greenscreen = EMF2_peaks_summit %>% inner_join(EMF2_peaks_greenscreen %>% dplyr::select(name))

# Extend the peaks to +/- 300bp
EMF2_peaks_summit_greenscreen_extend = EMF2_peaks_summit_greenscreen %>% mutate(start_extend=start-300, end_extend=end+299) %>% dplyr::select(chr, start_extend, end_extend, name, value)

# Export as bed to be converted to FASTA (format should be chr01...) 
write.table(EMF2_peaks_summit_greenscreen_extend, file="data/ChIPseeker/EMF2_peaks_summit_greenscreen_extend.bed",row.names=FALSE,quote=FALSE,sep='\t')

# Import EMF2 and H3K27me3 overlapping peaks and modify their peak lenght as previously to reach 600bp
overlap_genes_EMF2_complete = read_delim('data/ChIPseeker/overlap_genes_EMF2_complete.bed') # import
EMF2_peaks_summit_overlap_genes_complete = EMF2_peaks_summit %>% inner_join(overlap_genes_EMF2_complete %>% dplyr::select(name)) # combine
EMF2_peaks_summit_overlap_genes_complete_extend = EMF2_peaks_summit_overlap_genes_complete %>% mutate(start_extend=start-300, end_extend=end+299) %>% dplyr::select(chr, start_extend, end_extend, name, value) # extend

# Export as bed to be converted to FASTA (format should be chr01...) 
write.table(EMF2_peaks_summit_overlap_genes_complete_extend, file="data/ChIPseeker/EMF2_peaks_summit_overlap_genes_complete_extend.bed",row.names=FALSE,quote=FALSE,sep='\t')
```
Here is two files motif discovery can be assessed:
- EMF2_peaks_summit_greenscreen_extend.bed = All the EMF2 greenscreen peaks with a 600bp length
- EMF2_peaks_summit_overlap_genes_complete_extend.bed = All the EMF2 greenscreen peaks overlapping with H3K27me3 with a 600bp length
Now let's look at the binding profile of gene locus and create a control random genomic regions that have same lenght and same distribution over genes (as in [Winter et al 2011 Dev cell](XXX))
```R
library(ChIPseeker)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Distribution of EMF2_peaks_summit_overlap_genes_complete_extend
txdb <- makeTxDbFromBiomart(biomart="plants_mart",
                            dataset="osativa_eg_gene",
                            host="plants.ensembl.org")
                           

# Convert chr name as numerical value as in the txdb file

chr_label <- data.frame (Chr  = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"),
                  chr = c(1:12)
                  )

EMF2_peaks_summit_overlap_genes_complete_extend_ChrNumeric = EMF2_peaks_summit_overlap_genes_complete_extend %>% dplyr::rename(Chr=chr, start=start_extend, end=end_extend) %>% left_join(chr_label) %>% dplyr::select(chr, start, end, name, value)

# Create genomic range file
peaks.gr = makeGRangesFromDataFrame(EMF2_peaks_summit_overlap_genes_complete_extend_ChrNumeric, keep.extra.columns=TRUE)

# Annotate peak to genes
peakAnno  = annotatePeak(peaks.gr,tssRegion=c(-3000,500), TxDb=txdb)

# Check distribution plot relative to features
pdf('data/ChIPseeker/EMF2_peaks_summit_overlap_genes_complete_extend_distribution.pdf')
plotAnnoPie(peakAnno)
dev.off()

# Check distribution to gene
pdf('data/ChIPseeker/EMF2_peaks_summit_overlap_genes_complete_extend_distribution_gene.pdf')
plotPeakProf2(peak = peaks.gr, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 100,
              TxDb = txdb, ignore_strand = F)
dev.off()
```






Can try an aditional filter with motifs from "inducible genes" = very low express under the condition ChIP was done and higher express in some higher tissue = "inducible genes". --> Does TPM is ok?? 

Do not filter for promoter, keep them all!

For their length, keep summit +/- 300bp.
Generate a background genomic dataset not bound: random sequences same length and same distribution. Winter et al 2011 Dev. Cell LFY CHIP. 

 
 
 











see for next:
https://hbctraining.github.io/Intro-to-ChIPseq/lessons/12_functional_analysis.html
https://www.hdsu.org/chipatac2020/06_CHIP_PeakAnnotation.html



Follow [this](https://github.com/sklasfeld/ChIP_Annotation) method\
--> Proceed with the tutorial to understand what s going on






`








Minor issues:
- Mapping only work on node01 or node03
- 
