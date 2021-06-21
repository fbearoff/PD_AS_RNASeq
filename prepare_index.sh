#!/bin/bash

#https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

#supply from ensembl
genome="$HOME/hg19/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
transcriptome="$HOME/hg19/Homo_sapiens.GRCh38.cdna.all.fa.gz"

threads="$(grep -c ^processor /proc/cpuinfo)"
parentdir="$(dirname "$genome")"

grep "^>" <(gunzip -c $genome) | cut -d " " -f 1 > $parentdir/decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

cat $transcriptome $genome > $parentdir/gentrome.fa.gz

salmon index -t $parentdir/gentrome.fa.gz -d $parentdir/decoys.txt -p $threads -i $parentdir/salmon_index
