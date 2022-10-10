#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL

outdir="fastq/raw"
mkdir -p ${outdir}




x="EMF2_Rep1_1"
raw_f="SRR18596327_1.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi


x="EMF2_Rep1_2"
raw_f="SRR18596327_2.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi





x="EMF2_Rep2_1"
raw_f="SRR18596325_1.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi



x="EMF2_Rep2_2"
raw_f="SRR18596325_2.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi



x="IgG_1"
raw_f="SRR15663624_1.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi


x="IgG_2"
raw_f="SRR15663624_2.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi




x="H3K27me3_Rep1_1"
raw_f="SRR15663617_1.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="H3K27me3_Rep1_2"
raw_f="SRR15663617_2.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi



x="H3K27me3_Rep2_1"
raw_f="SRR15663616_1.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi


x="H3K27me3_Rep2_2"
raw_f="SRR15663616_2.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi







x="input_1"
raw_f="SRR15663629_1.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi



x="input_2"
raw_f="SRR15663629_2.fastq"
gzip_f="${outdir}/${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
        gzip -c ${raw_f} > ${gzip_f} \
                && rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
        echo "ERROR: Cannot Find File: ${raw_f}"
        exit
fi

