#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL

outdir="fastq/raw"
mkdir -p ${outdir}


x="SDG711RNAi_Rep1_2"
raw_f="SRR15663619_2.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi


x="SDG711RNAi_Rep2_1"
raw_f="SRR15663618_1.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi



x="SDG711RNAi_Rep2_2"
raw_f="SRR15663618_2.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="WT_Rep3_1"
raw_f="SRR15663623_1.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi


x="WT_Rep3_2"
raw_f="SRR15663623_2.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi


