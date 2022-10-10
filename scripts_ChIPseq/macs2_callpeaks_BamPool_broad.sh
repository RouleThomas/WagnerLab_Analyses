#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=20G
#SBATCH --nodelist=node01

raw_bam_dir="mapped/chip"
downsamp_bam_dir="${raw_bam_dir}/downsample"
bam_suffix="dupmark.sorted.bam"

macs2_out="data/macs2_out/chipPeaks"
gs_regions="../GreenScreen/rice/GreenscreenProject/data/macs2_out/inputControls/qval10/gs_merge10000bp_call10_20inputs.bed" #10kb merge greenscreen



# average basepair q-value threshold (log10)
q=10

# make macs2 output directory
mkdir -p ${macs2_out}/broad_noMask_qval${q}
mkdir -p ${macs2_out}/broad_gsMask_qval${q}





# Call peaks on both EMF2 replicates and filtered out overlaping peaks with greenscreen
macs2 callpeak -t mapped/chip/downsample/EMF2_Rep1.dupmark.sorted.bam mapped/chip/downsample/EMF2_Rep2.dupmark.sorted.bam \
  -c mapped/chip/downsample/input.dupmark.sorted.bam  \
  -f BAMPE --keep-dup auto \
  --nomodel -g 373128865 \
  --outdir ${macs2_out} -n EMF2_pool --broad
                                                                                                 
awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/EMF2_pool_peaks.broadPeak > ${macs2_out}/broad_noMask_qval${q}/EMF2_pool_peaks.broadPeak                                             
    
bedtools intersect -v -wa \
  -a ${macs2_out}/broad_noMask_qval${q}/EMF2_pool_peaks.broadPeak \
  -b ${gs_regions} > ${macs2_out}/broad_gsMask_qval${q}/EMF2_pool_peaks.broadPeak                                             
                                                 
# Call peaks on both H3K27me3 replicates and filtered out overlaping peaks with greenscreen
macs2 callpeak -t mapped/chip/downsample/H3K27me3_Rep1.dupmark.sorted.bam mapped/chip/downsample/H3K27me3_Rep2.dupmark.sorted.bam \
  -c mapped/chip/downsample/input.dupmark.sorted.bam  \
  -f BAMPE --keep-dup auto \
  --nomodel -g 373128865 \
  --outdir ${macs2_out} -n H3K27me3_pool --broad
                                                                                                 
awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/H3K27me3_pool_peaks.broadPeak > ${macs2_out}/broad_noMask_qval${q}/H3K27me3_pool_peaks.broadPeak                                             
    
bedtools intersect -v -wa \
  -a ${macs2_out}/broad_noMask_qval${q}/H3K27me3_pool_peaks.broadPeak \
  -b ${gs_regions} > ${macs2_out}/broad_gsMask_qval${q}/H3K27me3_pool_peaks.broadPeak    


