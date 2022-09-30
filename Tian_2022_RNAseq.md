# Re-analyses of the Tian et al 2022 RNAseq using:

- trimming: XXX
- mapping: XXX
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
fasterq-dump SRR18596327 -S

```
## 2. Start analyses ##

