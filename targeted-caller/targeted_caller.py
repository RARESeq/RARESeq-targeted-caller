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
parser.add_argument('-fusiontarget', metavar='', type=str, help='<file> List of target fusions. A tab-delimited file containing on the first 2 columns gene symbols for the fusion partners.')
parser.add_argument('-fusionoutput', metavar='', type=str, help='<file> STAR fusion output. Should be a file called "star-fusion.fusion_candidates.preliminary.wSpliceInfo.wAnnot".')
parser.add_argument('-background', metavar='', type=str, help='<file> Pickle file containg the background stats, generated using targeted-background.py script. ')
parser.add_argument('-g', metavar='', type=str, help='<file> Genome file used to map the sample. fai file must exist. (' + str(parser.get_default('g')) + ')')
parser.add_argument('-d', metavar='', type=int, help='<int> Minimum depth to consider position viable. (' + str(parser.get_default('d')) + ')')
parser.add_argument('-bt', metavar='', type=str, help='<str> The bedtools command.')
args = parser.parse_args() 

#Save input arguments
output_dir = args.output
depth = args.d
genome = args.g
target = args.target
background_file = args.background 
freqs = args.freqs
fusiontarget = args.fusiontarget
fusionoutput = args.fusionoutput

if output_dir == None:    
    output_dir = os.path.dirname(freqs[0])
    print("Setting output dir to default: " + output_dir)
os.makedirs(output_dir, exist_ok = True)

def calculate_posterior(nr, depth, alpha, beta):
    p = 1 - betabinom.cdf(nr - 1, depth, alpha, beta)
    return p

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


stats, genes = pickle.load(open(background_file, "rb"))
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
   

g = open(os.path.join(output_dir, "calls_SNV_idels.txt"), "w")
g.write("CHR\tPOS\tREF\tVAR\tDEPTH\tNR\tGENE\tPROB\tP_NR\tTYPE\tBKG_MAX_AF\n")

#Parse line by line on freq files
for freqfile in [x for x in freqs if '.freq.' in x]:
    print(f"Calling mutations in file: {freqfile}")        
    for line in open(freqfile,'r'):
        tokens = line.strip().split('\t')
        #Skip header
        if tokens[0] == 'CHR': continue
        #Skip low depth and snp
        if int(tokens[2]) < depth: continue

        pid = tokens[0]+":"+tokens[1]+":"+tokens[3]+":"
        
        if not pid in target_dict:     
            continue

        #Parse each type
        for bp,col in zip('ACTG',[6,8,10,12]):
            if not pid + bp in target_dict[pid]:
                continue            
            nr = int(tokens[col])+int(tokens[col+1])
            d = int(tokens[2])
            #afs = nr/d
            bkg_nr = stats[pid][pid+bp]["NR"]
            bkg_d = stats[pid][pid+bp]["DEPTH"]
            bkg_af = stats[pid][pid+bp]["AF"]
            alpha = stats[pid][pid+bp]["ALPHA"]
            beta = stats[pid][pid+bp]["BETA"]
            alpha_init = stats[pid][pid+bp]["ALPHA_INIT"]
            beta_init = stats[pid][pid+bp]["BETA_INIT"]
            
            gene_name = genes[pid + bp]
                       
            if alpha_init > 0:
                p = calculate_posterior(nr, d, alpha_init, beta_init)
                prob = alpha_init / (alpha_init + beta_init)
                mx = max(bkg_af)
                #print(alpha_init, beta_init, alpha, beta, p, nr, d)
                #if p < 0.05:                
                #    print(alpha_init, beta_init, alpha, beta, p, nr, d, gene_name, freqfile)
            else:
                prob = 0
                p = 0
                mx = 0         
            #print(("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(tokens[0], tokens[1], tokens[3], bp, tokens[2], nr, gene_name, p, p_af, z, mn, sd, str(len(bkg)))))
            g.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format(tokens[0], tokens[1], tokens[3], bp, tokens[2], nr, gene_name, prob, p, "SNV", mx))
        

#Parse line by line in indels files
for indelfile in [x for x in freqs if '.indels.' in x]:
    #For each line
    print(f"Calling indels in file: {indelfile}")        
    for line in open(indelfile,'r'):       
        tokens = line.strip().split('\t')
        #Skip header
        if tokens[0] == 'CHR': continue 
        #Skip low depth and snp
        if len(tokens) < 2: continue
        if int(tokens[2]) < depth: continue
                
        pid = tokens[0]+":"+tokens[1]+":"+tokens[3]+":"        
        short_pid = tokens[0]+":"+tokens[1]
        #print(line)
        type = "NA"
        if short_pid in special_target_dict:        
            #print(special_target_dict[short_pid])
            n = len(tokens[4]) - 1
            key1  = tokens[4][0] + "*"
            key2  = tokens[4][0] + str(n)
            if key1 in special_target_dict[short_pid]:
                #print("Variant called by generic rule 1. Variant info:\n" + line)
                type = "Indel"
                gene_name = genes[short_pid + ":" + key1]
            else:
                if key2 in special_target_dict[short_pid]:
                    #print("Variant called by generic rule 2. Variant info:\n" + line)
                    type = "Indel"
                    gene_name = genes[short_pid + ":" + key2]
            
        
        if  pid in target_dict and pid + tokens[4] in target_dict[pid]:                     
            gene_name = genes[pid + bp]
            type = "Indel"
        
        if type == "NA":                
            continue

        nr = int(tokens[5]) + int(tokens[6])
        afs = nr/int(tokens[2])
        d = int(tokens[2])
        bp = tokens[4] 
        p = 0
        prob = 1
        
        g.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format(tokens[0], tokens[1], tokens[3], bp, tokens[2], nr, gene_name, prob, p, type, "NaN")) 
    
g.close()

if fusiontarget == None or fusionoutput == None or not os.path.exists(fusiontarget) or not os.path.exists(fusionoutput):
    print("Fusion target or STAR fusion output not provided. Not calling fusions. Fusion parameters provided: -fusiontarget=", fusiontarget, " -fusionoutput", fusionoutput)
else:
    ground_truth_fusions = {}
    for line in open(fusiontarget, "r").readlines():
        tokens = line.strip().split("\t")
        if tokens[0]=="GENE1":
            continue
        
        if tokens[0] not in ground_truth_fusions:
            ground_truth_fusions[tokens[0]] = set()
        ground_truth_fusions[tokens[0]].add(tokens[1])
        if tokens[1] not in ground_truth_fusions:                    
            ground_truth_fusions[tokens[1]] = set()
        ground_truth_fusions[tokens[1]].add(tokens[0])

    #print(ground_truth_fusions)
    g = open(os.path.join(output_dir, "calls_fusions.txt"), "w")
    
    for line in open(fusionoutput, "r").readlines():        
        tokens = line.strip().split("\t")
        if tokens[0] == "#FusionName":
            g.write("\t".join(tokens) + "\n")
            continue

        partners = tokens[0].split("--")    
        #print(partners)     
        if partners[0] in ground_truth_fusions and (partners[1] in ground_truth_fusions[partners[0]] or "?" in ground_truth_fusions[partners[0]]) or \
           partners[1] in ground_truth_fusions and (partners[0] in ground_truth_fusions[partners[1]] or "?" in ground_truth_fusions[partners[1]]):
            #print(line)
            g.write("\t".join(tokens) + "\n")
    g.close()
