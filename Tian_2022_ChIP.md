# Re-analyses of the Tian et al 2022 ChIP using the greenscreen pipeline
- Tian et al 2022 paper [here]([https://pages.github.com/](https://academic.oup.com/plcell/article/34/8/2969/6580212#367433007)
- Klasfeld et al 2022 paper [here](https://www.biorxiv.org/content/10.1101/2022.02.27.482177v1)

## 1. import raw fastq file ##
--> Setup prerequisets
```
srun --mem=20g --pty bash -l
cd Tian_2022TPC_ChIP
module load sratoolkit/2.11.2
```
--> Download files from the SRA (paired end mode):
*Files information can be found in ChIP analyses Tan et al 2022.xlsx#
```
fasterq-dump SRR18596327 -S
fasterq-dump SRR18596325 -S

```
