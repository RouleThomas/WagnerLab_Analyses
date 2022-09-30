# Re-analyses of the Tan et al 2022 RNAseq using:
[Paper](https://academic.oup.com/plcell/article-abstract/34/8/2969/6580212?redirectedFrom=fulltext&login=false#367433007)

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
fasterq-dump SRR15663633 -S
fasterq-dump SRR15663632 -S
fasterq-dump SRR15663619 -S
fasterq-dump SRR15663618 -S
fasterq-dump SRR15663623 -S
```
ongoing

## 2. Start analyses ##

