#!/usr/bin/python3

import os
import sys
import re
import argparse
import pickle
import numpy
import subprocess
from scipy.stats import betabinom
from scipy.optimize import minimize 

#Get input arguments from command line
parser = argparse.ArgumentParser(description='Get background SNVs to filter for variant caller. Will output all positions matching the given criteria in the list of freq files given.')
parser.set_defaults(d=0,g='/oak/stanford/groups/andrewg/users/bluca/RARE-Seq/pipeline/indexes/hg19_ERCC.fa',fusiontarget=None,fusionoutput=None,output = None)
parser.add_argument('-output', help='<file> Name of output directory.')
parser.add_argument('-freqs',metavar='', help='Any number of freq files to be read in unison.', nargs=argparse.REMAINDER)
parser.add_argument('-target', metavar='', type=str, help='<file> List of target point mutations and/or indels. A tab-delimited file containing on the first 5 columns: chromosome, position, reference, varaint, gene_name.')
parser.add_argument('-g', metavar='', type=str, help='<file> Genome file used to map the sample. fai file must exist. (' + str(parser.get_default('g')) + ')')
parser.add_argument('-bt', metavar='', type=str, help='<str> The bedtools command.')
args = parser.parse_args() 

#Save input arguments
output_dir = args.output
genome = args.g
target = args.target
freqs = args.freqs

os.makedirs(output_dir, exist_ok = True)

def run_command(command):
    try:
        # Execute the command without capturing stdout and stderr
        result = subprocess.run(command, shell=True)
        
        # Check if the command failed
        if result.returncode != 0:
            print(f"Command '{command}' failed with return code {result.returncode}")
            sys.exit(result.returncode)
    
    except FileNotFoundError:
        print(f"Command not found: {command}")
        sys.exit(1)
    except Exception as e:
        # Catch any other exceptions and print the error
        print(f"Unexpected error: {str(e)}")
        sys.exit(1)


genomeorder = os.path.join(os.path.dirname(output_dir), 'genome-order.sorting.txt')
run_command('cut -f 1-2 {0}.fai > {1}'.format(genome, genomeorder))
gchro = os.popen('cut -f 1 {0}.fai'.format(genome)).read().strip().split('\n')

sorted_target = target.replace(".txt", "").replace(".bed", "") + ".sorted.bed" 
if not os.path.exists(sorted_target):
    header = os.popen('head -1 ' + target).read().strip().split('\t')
    run_command('awk \'{{if($1 !~ "CHR")print $1 "\\t" $2 "\\t" $3 "\\t" $0;}}\' {0} | {1} sort -faidx {2}.fai -i - >{3}'.format(target, args.bt, genome, sorted_target))

alluchro = os.popen('cut -f 1 {0} | uniq'.format(sorted_target)).read().strip().split('\n')
uchro = []
for chro in gchro:
    if chro in alluchro:
        uchro.append(chro)

genes = {}
target_dict = {}  
for line in open(target, "r"):
    tokens = line.strip().split()
    if tokens[0] == "CHR":
        continue
    
    if tokens[0] + ":" + tokens[1] + ":" + tokens[3] + ":" in target_dict:
        target_dict[tokens[0] + ":" + tokens[1] + ":" + tokens[3] + ":" ].append(tokens[0] + ":" + tokens[1] + ":" + tokens[3] + ":" + tokens[4])
    else:
        target_dict[tokens[0] + ":" + tokens[1] + ":" + tokens[3] + ":" ] = [tokens[0] + ":" + tokens[1] + ":" + tokens[3] + ":" + tokens[4]]
    genes[tokens[0] + ":" + tokens[1] + ":" + tokens[3] + ":" + tokens[4]] = tokens[5]

special_target_dict = {}  
for line in open(target, "r"):
    tokens = line.strip().split()
    if tokens[0] == "CHR":
        continue
    
    if "*" in tokens[4] or re.match("[+-][0-9]+?", tokens[4]):
        for i in range(int(tokens[1]), int(tokens[2]) + 1):                    
            if tokens[0] + ":" + str(i) in special_target_dict:
                special_target_dict[tokens[0] + ":" + str(i)].append([tokens[4]])
            else:
                special_target_dict[tokens[0] + ":" + str(i)] = [tokens[4]]
            genes[tokens[0] + ":" + str(i) + ":" + tokens[4]] = tokens[5]
   

for freq in freqs:    
    print("Filtering input file: {0}".format(freq))
    filteredfreq = re.sub('\.txt$', '.filtered.txt', freq)    
    
    if ".freq." in filteredfreq:        
        os.system('sh targeted-caller/dump_reads_overlapping_variants.sh {0}'.format(os.path.dirname(filteredfreq))) 
    
    #print("Filtering freq file: ", freq)
    header = os.popen('head -1 ' + freq).read().strip().split('\t')
    open(filteredfreq, "w").write("\t".join(header) + "\n")
    
    for chro in uchro:
        #print("Filtering chromosome: " + chro)
        run_command('grep -E "^{0}\s" {1} | head -n 1 - >> {2}'.format(
                chro, freq, filteredfreq))                
        
        run_command('grep -E "^{0}\s" {1} | awk \'{{if(NR>1)print $1 "\\t" $2 "\\t" $2 "\\t" $0}}\' - | {2} intersect -a stdin -b {3} -sorted -g {4} | cut -f 4-{5} - | rev | uniq -f 1 | rev >> {6}'.format(
                chro, freq, args.bt, sorted_target, genomeorder, len(header)+3, filteredfreq))                

