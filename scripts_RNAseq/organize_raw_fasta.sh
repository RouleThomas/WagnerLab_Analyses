#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL

outdir="fastq/raw"
mkdir -p ${outdir}



x="WT_Rep1_1"
raw_f="SRR15663633_1.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi


x="WT_Rep1_2"
raw_f="SRR15663633_2.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi


x="WT_Rep2_1"
raw_f="SRR15663632_1.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi


x="WT_Rep2_2"
raw_f="SRR15663632_2.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi



x="SDG711RNAi_Rep1_1"
raw_f="SRR15663619_1.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi



