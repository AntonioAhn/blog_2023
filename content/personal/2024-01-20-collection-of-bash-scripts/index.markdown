---
title: collection of bash scripts
author: ''
date: '2024-01-20'
slug: []
categories: []
tags: []
---


Just a collection of miscellaneous bash scripts


# macs2 with in input from 4 batches of bam files with 4 different control inputs


```bash
#!/bin/bash
#SBATCH --job-name=MACS2-job
#SBATCH -N 1
#SBATCH --mem=12G
#SBATCH --cpus-per-task=12
#SBATCH -t 12:00:00
#SBATCH --mail-user=antonio.ahn@petermac.org
#SBATCH --mail-type=end
#SBATCH --output=MACS-%j.out
#SBATCH --error=MACS-%j.err

# saved to macs2_withinput.sh /path/scripts

# load MACS
module load macs/2.2.7.1
module load bedtools/2.27.1

# path variable
# move to folder with BAM files
input_dir="/path/results/bam_files/final_bams/hg38"
output_files="/path/results/macs2/withinput/macs2_q0.05"
cd $input_dir

# make directory of macs2
mkdir -p $output_files


# list the 5 batches and the 5 control IgG's
# The printf statement will prepend the value with leading zeros so the width is 3, for example ("4" becomes "004").
batch_1=($(printf "%03d-MCF7-*.bam " {1..17}))
batch_2=($(printf "%03d-MCF7-*.bam " {19..43}))
batch_3=($(printf "%03d-MCF7-*.bam " {45..62}))
batch_4=($(printf "%03d-MCF7-*.bam " {64..87}))
# batch_5=($(printf "%03d-MCF7-*.bam " {89..105}))
controls=("018-MCF7-Par-DMSO-R1-Rabbit-IgG-native_S9_hg38_noecoli_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam" "044-MCF7-Par-DMSO-R1-Rabbit-IgG-native_S38_hg38_noecoli_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam" "062-MCF7-Par-DMSO-R2-Rabbit-IgG-native_S58_hg38_noecoli_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam" "063-MCF7-Par-DMSO-R1-Mouse-IgG2a-native_S59_hg38_noecoli_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam" "088-MCF7-Par-DMSO-R1-Rabbit-IgG-mod-fix_S86_hg38_noecoli_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam")


for i in {1..4};
do
# Determine the current batch array
current_batch="batch_$i[@]"

# Echo the batch number
echo "Processing Batch number $i"

# Loop through the current batch
for file in "${!current_batch}"; do

base=`basename ${file} _chrmrm_blrm_duprm_unmappedrm_multimaprm.bam`
echo "sample name is $base"
echo "control name is ${controls[$((i-1))]}"

macs2 callpeak -t $file \
-c ${controls[$((i-1))]} \
-f BAMPE -g hs \
-q 0.05 \
-n $base \
--outdir $output_files;

done;
done
```



# cat pairs of bam files


```bash
filepairs=($(ls *R1*.fastq.gz) $(ls *R2*.fastq.gz))

# this shows that the order is correct. tbh im not sure if i need to do the basenames here? 
# also it doesnt matter if run1 shows first then run2 (or the other way round) as long as they are the correct pairs. THey will be cat'ed
for (( idx=0 ; idx<${#filepairs[@]} ; idx+=2));
do
replicate_1=${filepairs[idx]}
replicate_2=${filepairs[idx+1]}

base_1=`basename ${filepairs[idx]} _001_run1.fastq.gz`
base_2=`basename ${filepairs[idx+1]} _001_run2.fastq.gz`

outputname=$(echo $replicate_1 | sed -e 's/001_run1/001_cat/g' -e 's/001_run2/001_cat/g' )

echo $replicate_1
echo $replicate_2
#echo $base_1
#echo $base_2
echo $outputname;

done

# now cat the files 
for (( idx=0 ; idx<${#filepairs[@]} ; idx+=2));
do
replicate_1=${filepairs[idx]}
replicate_2=${filepairs[idx+1]};

outputname=$(echo $replicate_1 | sed -e 's/001_run2/001_cat/g' -e 's/001_run1/001_cat/g' )

cat $replicate_1 $replicate_2 > ${outputname};

done

```



# deep tools heatmap, color up and dwon separately


```bash

#!/bin/bash
#SBATCH --job-name=deeptools-job
#SBATCH -N 1
#SBATCH --partition=prod_med
#SBATCH --mem=8G
#SBATCH --cpus-per-task=15
#SBATCH -t 10:00:00
#SBATCH --mail-user=antonio.ahn@petermac.org
#SBATCH --mail-type=end
#SBATCH --output=DT_heatmap-%j.out
#SBATCH --error=DT_heatmap-%j.err

# script to make one compute matrix files for separate bigwig files (DMSO and LY treated in this case)
# saved as heatmap.sh in /path/R_analysis/4.diffbind_overcome_resistance/1.overcome_resist_LYvsDMSO/2.annotations/deeptools_heatmap/scripts

module load deeptools/3.5.0

BW_dir="/path/results/bigwig"
BED_dir="/path/R_analysis/4.diffbind_overcome_resistance/1.overcome_resist_LYvsDMSO/2.annotations/data/bed"
output_dir="/path/R_analysis/4.diffbind_overcome_resistance/1.overcome_resist_LYvsDMSO/2.annotations/deeptools_heatmap/output"

mkdir -p $output_dir
set -xe

for i in up down;
do

if [[ $i == "up" ]]; then
 color=Reds
 echo $color
else
 color=Blues
 echo $color
fi

computeMatrix reference-point --referencePoint center \
-b 1000 -a 1000 \
-R $BED_dir/ParLY_vs_ParDMSO_promoter_${i}.bed $BED_dir/ParLY_vs_ParDMSO_nonpromoter_${i}.bed \
-S $BW_dir/Sample-1-DMSO-Rep-1_S10.bw \
$BW_dir/Sample-9-Cryo-Par-DMSO-Rep-2_S14.bw \
$BW_dir/Sample-2-Par-LY-Rep-1_S2.bw \
$BW_dir/Sample-10-Par-LY-Rep-2_S10.bw \
$BW_dir/Sample-5-LYR-DMSO-Rep-1_S5.bw \
$BW_dir/Sample-13-Cryo-LYR-DMSO-Rep-2_S13.bw \
$BW_dir/Sample-6-LYR-ARC-Rep-1_S12.bw \
$BW_dir/Sample-14-LYR-ARC-Rep-2_S14.bw \
$BW_dir/Sample-7-LYFR-DMSO-Rep-1_S7.bw \
$BW_dir/Sample-15-Cryo-LYFR-DMSO-Rep-2_S15.bw \
$BW_dir/Sample-8-LYFR-ARC-Rep-1-8-1_S13.bw \
$BW_dir/Sample-16-LYFR-ARC-Rep-1-8-2_S16.bw \
--skipZeros \
--binSize 10 \
--missingDataAsZero \
--sortUsingSamples 3 \
--sortRegions descend \
-p 14 \
-o $output_dir/PAR_LYvsDMSO_peaks_overcomeresist_${i}.gz

plotHeatmap -m $output_dir/PAR_LYvsDMSO_peaks_overcomeresist_${i}.gz \
-out $output_dir/PAR_LYvsDMSO_peaks_overcomeresist_${i}.png  \
--colorMap $color \
--dpi 300 \
--sortUsingSamples 3 \
--refPointLabel center \
--samplesLabel PAR_DMSO_R1 PAR_DMSO_R2 PAR_LY_R1 PAR_LY_R2 \
LYR_DMSO_R1 LYR_DMSO_R2 LYR_CDK2i_R1 LYR_CDK2i_R2 \
LYFR_DMSO_R1 LYFR_DMSO_R2 LYFR_CDK2i_R1 LYFR_CDK2i_R2 \
--yAxisLabel "Accessibility" \
--regionsLabel promoter nonpromoter;

done

```



```bash
# different yMax depending on promoter or nonpromoter
for i in ParLY_vs_ParDMSO_*.bed;
do

if echo $i | grep -q '_nonpromoter_'; then
	number=60
  echo $number
else 
  number=100
  echo $number
fi

base=$(basename ${i} .bed | sed s/ParLY_vs_ParDMSO_//)

computeMatrix reference-point --referencePoint center \
-b 1000 -a 1000 \
-R $BED_dir/${i} \
-S $BW_dir/Sample-1-DMSO-Rep-1_S10.bw \
$BW_dir/Sample-9-Cryo-Par-DMSO-Rep-2_S14.bw \
$BW_dir/Sample-2-Par-LY-Rep-1_S2.bw \
$BW_dir/Sample-10-Par-LY-Rep-2_S10.bw \
--skipZeros \
--binSize 10 \
--missingDataAsZero \
--sortUsingSamples 3 \
--sortRegions descend \
-p 14 \
-o $output_dir/ParLY_vs_ParDMSO_${base}.gz

plotProfile -m $output_dir/ParLY_vs_ParDMSO_${base}.gz \
     -out $output_dir/ParLY_vs_ParDMSO_${base}.png \
     --plotHeight=15 \
     --plotWidth=15 \
     --plotType=lines \
     --perGroup \
     --colors blue blue red red \
     --yMax=$number \
     --yMin=0 \
     --samplesLabel PAR_DMSO_R1 PAR_DMSO_R2 PAR_LY_R1 PAR_LY_R2 \
     --yAxisLabel "Accessibility" \
     --plotTitle ${i}
```
