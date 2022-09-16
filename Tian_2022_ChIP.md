# Re-analyses of the Tian et al 2022 ChIP using the greenscreen pipeline
- Tian et al 2022 paper [here]([https://pages.github.com/](https://academic.oup.com/plcell/article/34/8/2969/6580212#367433007)
- Klasfeld et al 2022 paper [here](https://www.biorxiv.org/content/10.1101/2022.02.27.482177v1)

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
Submitted batch job 197661   \
Submitted batch job 197662   \



Minor issues:
- 
