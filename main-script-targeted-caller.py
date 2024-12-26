#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import shutil
sys.path.append("targeted-caller")
import demuxFunctions as demux
import time

parser = argparse.ArgumentParser(description='Targeted caller for the RARE-seq pipeline.')
parser.set_defaults(\
    genome='hg19_ERCC.fa',\
    background_samples="25",\
    threads=10,\
    seed=1234,\
    bedtools_command="bedtools")

parser.add_argument('-i', '--input_sample_list', required=True, type=str,
                    help='<file> Path to a tab-delimited file containing the sample list to be analyzed. \
                    This file should have at least 2 columns (matching of column names is case-sensitive): \
                    Lane: the lane on which the sample was sequenced, \
                    ID: sample ID.')
parser.add_argument('-b', '--background_samples', required=True, type=str,
                    help='<int|file> How the control samples for background estimation are selected. \
                    It can take as input either an integer number indicating the number of unique cotrol donors\
                    to be randomly selected from file provided in argument -i (sample IDs matching ^CTR*),\
                    or a file containing a subset of the rows provided in argument -i.')
parser.add_argument('-p', '--parsed_input', required=True, type=str, help='<file> Path containing the parsed lane outputs.')
parser.add_argument('-wm', '--whitelist_mutation_indel', required=True, type=str,
                    help='<file> Path to a tab-delimited whitelist file for point mutations and indels. \
                    This file should contain at least the following columns: \
                    CHR: chromosome, \
                    START: alteration start, \
                    END: alteration end (same as start for point mutations), \
                    REF: the reference allele, \
                    VAR: the variant allele. In case of indels, it should contain "-*" or "+*", \
                          for deletions and insertions respectively, \
                    GENE: The gene in which the alteration occurs.')
parser.add_argument('-wf', '--whitelist_fusion', required=True, type=str,
                    help='<file> Path to a tab-delimited whitelist file for gene fusions. \
                    This file should contain at least two columns, the gene symbols of fusion partners.')
parser.add_argument('-o', '--output_path', required=True, type=str,
                    help='<file> Output path for the current run.')
parser.add_argument('-g', "--genome", type=str, \
                    help='<file> Genome file used to map the sample. fai file must exist. (' + str(parser.get_default('g')) + ')')
parser.add_argument('-t', '--threads', required=False, type=int, 
                    help='<int> Number of threads.')
parser.add_argument('-s', '--seed', required=False, type=int, 
                    help='<int> The random seed.')
parser.add_argument('-bt', '--bedtools_command', required=False, type=str, 
                    help='<str> The bedtools command.')

args = parser.parse_args()
print(args.input_sample_list)

#Argument validation
if not os.path.isfile(args.input_sample_list):
    sys.exit(f"Error: The sample list file '{args.input_sample_list}' does not exist.")

if not os.path.isfile(args.whitelist_mutation_indel):
    sys.exit(f"Error: The whitelist mutation/indel file '{args.whitelist_mutation_indel}' does not exist.")

if not os.path.isfile(args.whitelist_fusion):
    sys.exit(f"Error: The whitelist fusion file '{args.whitelist_fusion}' does not exist.")

if not os.path.isdir(args.output_path):
    os.makedirs(args.output_path, exist_ok=True)

#Making sure the ouput dirrectory exists
output_dir = os.path.join(args.output_path)
parsed_input_dir = args.parsed_input
os.makedirs(output_dir, exist_ok=True)
print("Saving the outputs in: " + output_dir)
print(args)

parsed_input = args.parsed_input
backgroud_input = os.path.join(args.parsed_input, "parsed_background_files")

try:      
    print ("(1) Generating background stats.") 
    freq_files = []  
    f = open(os.path.join(parsed_input_dir, "background_files", "training_samples.txt"), "r")
    for line in f.readlines():
        tokens = line.split("\t")            
        if tokens[4] != "Freq" and tokens[4] != '"Freq"':
            freq_files.append(tokens[4].replace('"', ''))
    arg = ['python3', 'targeted-caller/targeted-background.py', "-output", os.path.join(output_dir, "background.p"),\
        "-target", args.whitelist_mutation_indel, "-bt", args.bedtools_command, "-g", args.genome,  "-freqs"]        
    
    for f in freq_files:         
        fl =  os.path.join(parsed_input_dir, "background_files", os.path.basename(os.path.dirname(f))+".freq.paired.Q30.txt")        
        arg.append(fl)
    
    result = subprocess.run(arg, check=True)
except subprocess.CalledProcessError as e:
    sys.exit(f"An error occurred: {e}")

try:      
    print ("(2) Calling variants.") 
    cmds = []
    f = open(os.path.join(parsed_input_dir, "sample_info.txt"), "r")
    for line in f.readlines():
        tokens = line.split("\t")            
        if tokens[4] == "Freq" or tokens[4] == '"Freq"':        
            continue
        
        lane = tokens[0].replace('"', '')
        sample = tokens[1].replace('"', '')    
        freq_file = tokens[4].replace('"', '').replace('\.txt$', '.filtered.txt')    
        indel_file = tokens[5].replace('"', '').replace('\.txt$', '.filtered.txt')    
        fusion_file = tokens[6].replace('"', '')

        os.makedirs(f"{output_dir}/call_data/{lane}_{sample}", exist_ok = True)
        if os.path.exists(fusion_file):        
            cmd = f"python3 targeted-caller/targeted_caller.py -output {output_dir}/call_data/{lane}_{sample} -background {os.path.join(output_dir, 'background.p')} \
            -bt {args.bedtools_command} -g {args.genome}\
            -target {args.whitelist_mutation_indel} -fusiontarget {args.whitelist_fusion} -fusionoutput {fusion_file} -freqs {freq_file} {indel_file}"
        else:        
            cmd = f"python3 targeted-caller/targeted_caller.py -output {output_dir}/call_data/{lane}_{sample} -background {os.path.join(output_dir, 'background.p')} \
            -bt {args.bedtools_command} -g {args.genome}\
            -target {args.whitelist_mutation_indel} -fusiontarget {args.whitelist_fusion} -freqs {freq_file} {indel_file} "

        cmds.append([cmd.replace("python3 ", "")])            
        
    demux.forkpoollang("python3", cmds, args.threads)
except subprocess.CalledProcessError as e:
    sys.exit(f"An error occurred: {e}")

#Aggregating results
try:      
    print ("(3) Aggregating the results.")     
    subprocess.run(['Rscript', 'targeted-caller/aggregate_results.R', parsed_input_dir, output_dir], check = True)    
except subprocess.CalledProcessError as e:
    sys.exit(f"(5) An error occurred: {e}")

#Plotting results
try:      
    print ("(4) Plotting the results.")     
    subprocess.run(['Rscript', 'targeted-caller/plot_results.R', parsed_input_dir,  output_dir], check = True)    
except subprocess.CalledProcessError as e:
    sys.exit(f"(6) An error occurred: {e}")

