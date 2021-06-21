#!/bin/bash

################################################################################
#Created by Frank Bearoff                                                      #
#Drexel University College of Medicine Department of Microbiology & Immunology #
################################################################################

read_files="/home/frank/data/FB38_June2021_30-506016007/00_fastq"
destination_directory="/home/frank/R_projects/FB38_June2021_30-506016007/salmon_output"
index="/home/frank/hg19/salmon_index"

#quantify transcripts against index
for sample in $read_files/*_R1_001.fastq.gz; do
    base_name="${sample##*/}"
    echo "Processing sample ${base_name%%_R1_001*}"
    salmon quant \
        -i "$index" \
        --gcBias \
        -l A \
        -1"${sample}" \
        -2 "${sample%%_R1*}_R2_001.fastq.gz" \
        -o "$destination_directory"/quants/"${base_name%%_R1*}"
    mapped=$(grep "percent" "$destination_directory"/quants/"${base_name%%_R1*}"/aux_info/meta_info.json|cut -d : -f 2|cut -d , -f 1)
    echo ${base_name%%_R1*} $mapped >> $destination_directory/mapped_percent.txt
done

