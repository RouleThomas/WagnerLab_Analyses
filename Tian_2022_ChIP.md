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







`








Minor issues:
- Mapping only work on node01 or node03
