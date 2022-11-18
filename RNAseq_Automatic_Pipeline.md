**Objective:** Design a python script that take as input raw fastq/download fastq and output DGE with DESeq2\
Let's modify the [Greenscreen demo script](https://github.com/sklasfeld/GreenscreenProject/blob/main/scripts/greenscreenPipeline.py) for this purpose.\

Let's anaylze several RNAseq datasets from [Zhang et al 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2860166/); let's take only the seedling shoot and root data to train:

GSM417534 	seedling shoot RNA-Seq data SRX017635
GSM417535 	seedling root RNA-Seq data  SRX017636
GSM417536 	tillering leaf RNA-Seq data SRX017637
GSM417537 	flowering panicle RNA-Seq data  SRX017638
GSM417538 	flowering leaf RNA-Seq data SRX017639
GSM417539 	filling panicle RNA-Seq data  SRX017640

