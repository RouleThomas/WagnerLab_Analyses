# Re-analyses of the Tian et al 2022 ChIP using the greenscreen pipeline
- Tian et al 2022 paper [here](https://academic.oup.com/plcell/article/34/8/2969/6580212#367433007)
- Klasfeld et al 2022 paper [here](https://www.biorxiv.org/content/10.1101/2022.02.27.482177v1)\
Tian et al used IgG as control for EMF2b ChIP, let's try using input instead

## 1. import raw fastq file ##
--> Setup prerequisets and Greenscreen conda environment
```
srun --mem=20g --pty bash -l
cd Tian_2022TPC_ChIP
module load sratoolkit/2.11.2
module load Anaconda/2019.10
conda activate CondaGS
```
Installed sound effect to be notify on command runnin more than 5seconds by following [this](https://github.com/c0rp-aubakirov/notify_after_command_executed/)
```
 git clone https://github.com/c0rp-aubakirov/notify_after_command_executed.git
 cd notify_after_command_executed/
 echo "source $(pwd)/postexec_notify" >> ~/.bashrc
 ```
 Restart terminal and tested if notification appear after 5sec ```sleep 6``` : \
 @FAIL:  Does not work when I restart, mention "postexec_notify: No such file or directory"
 
--> Download files from the SRA (paired end mode):\
*Files information can be found in ChIP analyses Tan et al 2022.xlsx*
```
fasterq-dump SRR18596327 -S
fasterq-dump SRR18596325 -S
fasterq-dump SRR15663624 -S
fasterq-dump SRR15663617 -S
fasterq-dump SRR15663616 -S
fasterq-dump SRR15663629 -S
```
## 2. Start analyses ##
*20220915*\
--> Compress/tidy/rename file with .sh script (script edited with ```nano```)
```
sbatch scripts/organize_raw_fasta.sh
```
Submitted batch job 197628=DONE\
*20220916*\
--> FASTQC to check adaptors type
```
sbatch scripts/fastqc_raw.sh
```
Submitted batch job 197643=DONE\
--> Trimming (some have adapters; some not; see *Tan et al 2022.xlsx*)
```
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
```
sbatch scripts/mapped_1.sh
```
Submitted batch job 197710 (FIE)=DONE\
Submitted batch job 197712 (H3K27me3)=DONE\
Submitted batch job 197714 (IGG)=DONE\
Submitted batch job 197715 (input)=DONE\
Start/End 1.20pm\10pm End \ 

--> Downsample replicate (for H3K27me3 and EMF2 2 Replicates)
```
sbatch scripts/downsampleBam.sh
```
Submitted batch job 197750\
Start/End 11.48am/=FAIL\
**output:**\
*[SEVERE][Biostar145820]Problem writing temporary file file:///tmp/sortingcollection.3596660749531177790.tmp.  Try setting TMP_DIR to a file system with lots of space.
htsjdk.samtools.util.RuntimeIOException: Problem writing temporary file file:///tmp/sortingcollection.3596660749531177790.tmp.  Try setting TMP_DIR to a file system with lots of space.*\
**troubleshoots:**\
- Try create a tmp folder in working directory to store temporary files and edit java comand (follow [this](https://www.biostars.org/p/42613/)) 
```
mkdir tmp
java -Djava.io.tmpdir=`pwd`/tmp TMP_DIR=`pwd`/tmp -jar COMAND
```
Submitted batch job 197762=FAIL\
**output:**\
*Error: Could not find or load main class TMP_DIR=.home.roule.Tian_2022TPC_ChIP.tmp*\
**troubleshoots:**\
- Try direct to tmp folder directly (follow [this](https://www.biostars.org/p/42613/)) 
```
java -Djava.io.tmpdir=tmp TMP_DIR=tmp -jar COMAND
```
Submitted batch job 197763=FAIL\
**output:**\
Error: Could not find or load main class TMP_DIR=tmp
**troubleshoots:**\
- Try to put the tmp-related argument at the end of the command
Submitted batch job 197764=WORK BUT NEW FAIL\
**output:**\
*[SEVERE][Biostar145820]There was an error in the command line arguments. Please check your parameters : Expected one or zero argument but got 3 : [mapped/chip/EMF2_Rep2.dupmark.sorted.bam, -Djava.io.tmpdir=tmp, TMP_DIR=tmp]
com.github.lindenb.jvarkit.lang.JvarkitException$CommandLineError: There was an error in the command line arguments. Please check your parameters : Expected one or zero argument but got 3*
**troubleshoots:**\
- Try to re-order correctly the different arguments and add ```-Xmx2g``` to allow 2G max of temporary files
```
 java -Xmx2g -Djava.io.tmpdir=tmp -jar ${JVARKIT_PATH}/dist/biostar145820.jar \
     --seed ${seed} -n ${min_val} \
     -o ${out_dir}/${samp}_Rep${rep}.dupmark.bam \
     ${in_dir}/${samp}_Rep${rep}.dupmark.sorted.bam \
     TMP_DIR=tmp
```
Submitted batch job 197767=FAIL\
**output:**\
*[SEVERE][Biostar145820]There was an error in the command line arguments. Please check your parameters : Expected one or zero argument but got 2 : [mapped/chip/EMF2_Rep2.dupmark.sorted.bam, TMP_DIR=tmp]\*
**troubleshoots:**\
- Remove the ```TMP_DIR=tmp```
Submitted batch job 197769=DONE\
Last=2hrs for 2 sample with 2 replicates\
**Troubleshootings Conclusion:** The java syntax was not good; For me ```-Xmx2g -Djava.io.tmpdir=tmp``` (store temporary file to tmp folder and allow 2G max) needs to be in between ```java``` (call java) and ```-jar``` (JAR file is used: biostar145820.jar); Then other arguments ```-n``` minimum number of read in one of the replicate; ```--seed 42``` (because that is the number of the universe, lol, could be any number... That is to generate randomness and reproducibility as a subset of the reads are taken)

--> Call peak IP to input control\
```
module load Anaconda/2019.10
conda activate CondaGS
sbatch scripts/macs2_callpeaks.sh
```
Submitted batch job 197779=FAIL\
**output:**\
*INFO  @ Tue, 20 Sep 2022 17:43:30: #3 Call peaks for each chromosome... 
ValueError: cannot resize this array: it does not own its data*
**troubleshoots:**\
Seems OK to ignore the error. Nevertheless, seems to be cause by not enough space in temporary shared folder, let's try adding one with the option ```--tempdir PATH```
```
module load Anaconda/2019.10
conda activate CondaGS
sbatch scripts/macs2_callpeaks_test.sh
```
Submitted batch job 197801=DONE but other FAIL\
No more error message but files are the same so that could have been ignore.\
**output:**\
*awk: `10' argument to `-v' not in 'var=value' form*
Need to add ```q=${q}``` after the -v argument, as follow: ```awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.narrowPeak > ${macs2_out}/noMask_qval${q}/${x}_peaks.narrowPeak``` --> Script XXX has been corrected





Minor issues:
- Mapping only work on node01 or node03
