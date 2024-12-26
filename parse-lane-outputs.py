#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import shutil
sys.path.append("targeted-caller")
import demuxFunctions as demux
import time

parser = argparse.ArgumentParser(description='Script that parses RARE-seq pipeline outputs into the format required by the targeted variat caller.')
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
parser.add_argument('-l', '--lane_path', required=True, type=str, help='<file> Path where the raw lane outputs are stored.')
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
os.makedirs(output_dir, exist_ok=True)
print("Saving the outputs in: " + output_dir)
print(args)

try:
    print ("(1) Collecting lane outputs")
    result = subprocess.run(['Rscript', 'targeted-caller/parse_lane_output.R', args.input_sample_list, args.lane_path, os.path.join(args.output_path, "parsed_lane_output")], check=True)    
except subprocess.CalledProcessError as e:
    sys.exit(f"An error occurred: {e}")

try:    
    if args.background_samples.isdigit():
        print ("(2) Selecting the background samples.")            
        result = subprocess.run(['Rscript', 'targeted-caller/select_training_samples_random.R', str(args.seed), str(args.background_samples), os.path.join(args.output_path, "parsed_lane_output"), os.path.join(args.output_path, "parsed_lane_output", "background_files")], check=True)
    else:
        print ("(2) Preparig the background samples.")            
        result = subprocess.run(['Rscript', 'targeted-caller/select_predefined_training_samples.R', str(args.seed), str(args.background_samples), os.path.join(args.output_path, "parsed_lane_output"), os.path.join(args.output_path, "parsed_lane_output", "background_files")], check=True)
except subprocess.CalledProcessError as e:
    sys.exit(f"An error occurred: {e}")

try:      
    print ("(3) Parsing background files.") 
    freq_files = []  
    f = open(os.path.join(output_dir, "parsed_lane_output", "background_files", "training_samples.txt"), "r")
    for line in f.readlines():
        tokens = line.split("\t")                 
        if tokens[4] != "Freq" and tokens[4] != '"Freq"':
            freq_files.append(tokens[4].replace('"', ''))    
    
    arg = ['python3', 'targeted-caller/targeted-background-filter-freq.py', "-output", os.path.join(output_dir, "parsed_lane_output", "background_files"),\
        "-target", args.whitelist_mutation_indel, "-bt", args.bedtools_command, "-g", args.genome, "-freqs"]        
    
    for f in freq_files:                        
        arg.append(f)
    
    result = subprocess.run(arg, check=True)
except subprocess.CalledProcessError as e:
    sys.exit(f"An error occurred: {e}")

try:   
    print ("(4) Parsing lane outputs.") 
    cmds = []
    f = open(os.path.join(args.output_path, "parsed_lane_output", "sample_info.txt"), "r")
    for line in f.readlines():
        tokens = line.split("\t")            
        if tokens[4] == "Freq" or tokens[4] == '"Freq"':        
            continue
        
        lane = tokens[0].replace('"', '')
        sample = tokens[1].replace('"', '')    
        freq_file = tokens[4].replace('"', '')
        indel_file = tokens[5].replace('"', '')
        fusion_file = tokens[6].replace('"', '')

        cmd = f"python3 targeted-caller/parse-input-files.py -output {output_dir}/parsed_lane_output \
        -bt {args.bedtools_command} -g {args.genome}\
        -target {args.whitelist_mutation_indel} -freqs {freq_file} {indel_file} "

        cmds.append([cmd.replace("python3 ", "")])            
    
    demux.forkpoollang("python3", cmds, args.threads)
except subprocess.CalledProcessError as e:
    sys.exit(f"An error occurred: {e}")
