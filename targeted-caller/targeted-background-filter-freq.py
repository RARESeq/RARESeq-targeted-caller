#!/usr/bin/python3

import os
import re
import argparse
import subprocess

# Get input arguments
parser = argparse.ArgumentParser(description='Filter background SNVs based on target mutations.')
parser.add_argument('-freqs', metavar='', nargs='+', required=True, 
                    help='Any number of freq files to be read in unison.')
parser.add_argument('-output_dir', metavar='', required=True, 
                    help='Output directory for filtered freq files.')
parser.add_argument('-target', metavar='', required=True, type=str, 
                    help='List of target point mutations and/or indels (tab-delimited).')
parser.add_argument('-g', metavar='', required=True, type=str, 
                    help='Genome file (with .fai index).')
parser.add_argument('-d', metavar='', type=int, default=0, 
                    help='Minimum depth to consider a position viable. Default: 0.')
parser.add_argument('-bt', metavar='', type=str, default='bedtools', 
                    help='Bedtools command (default: "bedtools").')
parser.add_argument('-force', action='store_true', 
                    help='Force re-filtering of freq files.')
args = parser.parse_args()

# Ensure output directory exists
filtered_freq_dir = os.path.join(args.output_dir)
os.makedirs(filtered_freq_dir, exist_ok=True)
genome = args.g
target = args.target

genomeorder = os.path.join(filtered_freq_dir, 'genome-order.sorting.txt')
subprocess.call('cut -f 1-2 {0}.fai > {1}'.format(genome, genomeorder), shell=True)
gchro = os.popen('cut -f 1 {0}.fai'.format(genome)).read().strip().split('\n')

sorted_target = target.replace(".txt", "").replace(".bed", "") + ".sorted.bed" 
if not os.path.exists(sorted_target):
    header = os.popen('head -1 ' + target).read().strip().split('\t')
    subprocess.call('awk \'{{if($1 !~ "CHR")print $1 "\\t" $2 "\\t" $3 "\\t" $0;}}\' {0} | {1} sort -faidx {2}.fai -i - >{3}'.format(target, args.bt, genome, sorted_target), shell=True)

alluchro = os.popen('cut -f 1 {0} | uniq'.format(sorted_target)).read().strip().split('\n')
uchro = []
for chro in gchro:
    if chro in alluchro:
        uchro.append(chro)

# Filter each frequency file
print("Filtering input files...")

for freq in args.freqs:
    output_subdir = os.path.join(filtered_freq_dir)
    os.makedirs(output_subdir, exist_ok=True)

    filteredfreq = os.path.join(output_subdir, os.path.basename(os.path.dirname(freq))+".freq.paired.Q30.txt")
    
    print(f"Filtering freq file: {freq}")
    header = os.popen('head -1 ' + freq).read().strip().split('\t')
    open(filteredfreq, "w").write("\t".join(header) + "\n")
    
    for chro in uchro:        
        subprocess.call('grep -E "^{0}\s" {1} | head -n 1 - >> {2}'.format(
                chro, freq, filteredfreq), shell=True)                        
        
        subprocess.call('grep -E "^{0}\s" {1} | awk \'{{if(NR>1)print $1 "\\t" $2 "\\t" $2 "\\t" $0}}\' - | {2} intersect -a stdin -b {3} -sorted -g {4} | cut -f 4-{5} - | rev | uniq -f 1 | rev >> {6}'.format(
                chro, freq, args.bt, sorted_target, genomeorder, len(header)+3, filteredfreq), shell=True)
    
print("Filtering complete.")
print("Filtered files saved to:", filtered_freq_dir)
