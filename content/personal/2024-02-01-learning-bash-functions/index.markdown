---
title: learning bash functions
author: ''
date: '2024-02-01'
slug: []
categories: []
tags: []
---


Currently learning bash functions. First script using it. 


```bash
#!/bin/bash
#SBATCH --job-name=deeptools-job
#SBATCH -N 1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=15
#SBATCH -t 10:00:00
#SBATCH --mail-user=antonio.ahn@petermac.org
#SBATCH --mail-type=end
#SBATCH --output=DT_heatmap-%j.out
#SBATCH --error=DT_heatmap-%j.err

# script to make one compute matrix files for separate bigwig files (DMSO and LY treated in this case)
# saved as deeptools_heatmap.sh in /scripts

module load deeptools/3.5.0

project_dir="/researchers/antonio.ahn"
BW_dir="${project_dir}/results/bigwig/BPM"
BED_dir="${project_dir}/R_analysis/2.peak_exploratory/data/bed"
output_dir="${project_dir}/R_analysis/2.peak_exploratory/deeptools_heatmap/output"

set -xe

# Function to process bigwig files
function process_files() {
    local file_array=("${!1}")  # Pass the array name as an argument
    local output_subdir="$2"
    local file_label="$2"
    
    mkdir -p $output_dir/$output_subdir

    computeMatrix reference-point --referencePoint center \
      -b 1000 -a 1000 \
      -R $BED_dir/consensus_regions.bed \
      -S "${file_array[@]}" \
      --skipZeros \
      --binSize 10 \
      --missingDataAsZero \
      --sortUsingSamples 3 \
      --sortRegions descend \
      -p 14 \
      --outFileSortedRegions $output_dir/$output_subdir/${file_label}_sorted.bed \
      -o $output_dir/$output_subdir/${file_label}_consensus_matrix.gz

    plotHeatmap -m $output_dir/$output_subdir/${file_label}_consensus_matrix.gz \
      -out $output_dir/$output_subdir/${file_label}.png  \
      --colorMap Reds \
      --dpi 300 \
      --sortUsingSamples 3 \
      --refPointLabel center \
      --samplesLabel Par_DMSO_1 Par_DMSO_2 Par_Abema_1 Par_Abema_2 sgRB1_DMSO_1 sgRB1_DMSO_2 sgRB1_Abema_1 sgRB1_Abema_2 \
      --yAxisLabel "${file_label} enrichment RPKM"
}

cd $BW_dir

# Process antibody_1_bw
antibody_1_bw=(*$(printf "%03d-MCF7-*bw " {45..52}))
process_files antibody_1_bw[@] "antibody_1"

# Process antibody_2_bw
antibody_2_bw=(*$(printf "%03d-MCF7-*bw " {53..60}))
process_files antibody_2_bw[@] "antibody_2"

# Process antibody_3_bw
antibody_3_bw=(*$(printf "%03d-MCF7-*bw " {19..26}))
process_files antibody_3_bw[@] "antibody_3"

# Process antibody_4_bw
antibody_4_bw=(*$(printf "%03d-MCF7-*bw " {1..8}))
process_files antibody_4_bw[@] "antibody_4"
```
