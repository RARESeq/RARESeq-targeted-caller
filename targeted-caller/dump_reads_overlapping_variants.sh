#!/bin/bash
output=$1

(samtools1.12 view -H $output/*.sorted.dualindex-deduped.sorted.bam; samtools1.12 view -L targeted-caller/MET14.bed $output/*.sorted.dualindex-deduped.sorted.bam) | samtools1.12 view -b -o $output/called_snv.cfRNA.bam -
samtools1.12 view -o $output/called_snv.cfRNA.sam $output/called_snv.cfRNA.bam
samtools1.12 index $output/called_snv.cfRNA.bam

cat $output/called_snv.cfRNA.sam | wc -l >$output/called_snv.cfRNA_METex14_skipping.sam 
cat $output/called_snv.cfRNA.sam | grep 3226N  - | wc -l >>$output/called_snv.cfRNA_METex14_skipping.sam 
cat $output/called_snv.cfRNA.sam | grep 3226N  - >$output/called_snv.cfRNA_METex14_skipping_raw.sam 
