#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL


macs2_out="data/macs2_out/chipPeaks"
gs_regions="../GreenScreen/rice/GreenscreenProject/data/macs2_out/inputControls/qval10/gs_merge10000bp_call10_20inputs.bed" #10kb merge greenscreen


# average basepair q-value threshold (log5)
q=10

# make macs2 output directory
mkdir -p ${macs2_out}/broad_noMask_qval${q}
mkdir -p ${macs2_out}/broad_gsMask_qval${q}


input_list=("IgG")

for x in "${input_list[@]}"; do
                               
awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.broadPeak > ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak                                             
    
bedtools intersect -v -wa \
  -a ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak \
  -b ${gs_regions} > ${macs2_out}/broad_gsMask_qval${q}/${x}_peaks.broadPeak                                             
                                                 
done                                                
