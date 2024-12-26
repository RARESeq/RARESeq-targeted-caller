#!/usr/bin/python3

import os
import sys
import re
import argparse
import pickle
import numpy
import subprocess
import random
from scipy.stats import betabinom
from scipy.optimize import minimize 

#Get input arguments from command line
parser = argparse.ArgumentParser(description='Get background SNVs to filter for variant caller. Will output all positions matching the given criteria in the list of freq files given.')
parser.set_defaults(d=0,f=False,g='hg19_ERCC.fa',\
    bt="bedtools")
parser.add_argument('-output', help='<file> Name of output pickle file. If not .p extention given one will be added')
parser.add_argument('-freqs',metavar='', help='Any number of freq files to be read in unison', nargs=argparse.REMAINDER)
parser.add_argument('-force',metavar='', type=bool, help='Force filtering of bam files.')
parser.add_argument('-target', metavar='', type=str, help='<file> List of target point mutations. A tab-delimited file containing on the first 6 columns: chromosome, start, end, reference, varaint, gene_name.')
parser.add_argument('-g', metavar='', type=str, help='<file> Genome file used to map the sample. fai file must exist. (' + str(parser.get_default('g')) + ')')
parser.add_argument('-d', metavar='', type=int, help='<int> Minimum depth to consider position viable. (' + str(parser.get_default('d')) + ')')
parser.add_argument('-bt', metavar='', type=str, required=False, help='<str> The bedtools command.')
args = parser.parse_args()

#Save input arguments
output = args.output
depth = args.d
genome = args.g
force = bool(args.force)
target = args.target
outfile = args.output


freqs = args.freqs
output_dir = os.path.dirname(output)

filtered_freq_dir = os.path.join(output_dir)
try:
    os.makedirs(filtered_freq_dir)
except:
    pass

def neg_log_likelihood(params, support_reads, total_reads):
    a, b = params
    nll = -numpy.sum(betabinom.logpmf(support_reads, total_reads, a, b))
    return nll

def fit_beta_binomial(counts, depths, alpha_init=1, beta_init=1):
    sorted_indices = numpy.lexsort((counts, depths))
    counts = numpy.array(counts)[sorted_indices]
    depths = numpy.array(depths)[sorted_indices]

    random.seed(1234)
    initial_params = [alpha_init, beta_init]
    bounds = [(0.01, None), (0.01, None)]
    result = minimize(neg_log_likelihood, initial_params, args=(counts, depths), bounds=bounds) 
    alpha, beta = result.x
    return alpha, beta

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

        
stats = {}
print("Print parsing freq files...")
#Parse line by line on freq files
for line in os.popen('python3 targeted-caller/read-freqs-by-line.py -p {0} '.format(depth) + ' '.join([x for x in freqs if '.freq.' in x])):
    #Get lines     
    lines = [x.split('\t') for x in line.split(':')]
    pid = lines[0][0]+":"+lines[0][1]+":"+lines[0][3]+":"

    if not pid in target_dict:     
        continue

    #Parse each type
    for bp,col in zip('ACTG',[6,8,10,12]):
        if not pid + bp in target_dict[pid]:
            continue

        afs = [(int(x[col])+int(x[col+1]))/int(x[2]) for x in lines]
        nr = [(int(x[col])+int(x[col+1])) for x in lines]
        d = [int(x[2]) for x in lines]   

        if sum(nr) > 0:
            alpha_init, beta_init = fit_beta_binomial(nr, d)                
            filtered_bkg_d = []
            filtered_bkg_nr = []
            for i in range(0, len(nr)):
                if nr[i] > 0:
                    filtered_bkg_nr.append(nr[i])
                    filtered_bkg_d.append(d[i])
            alpha_refined, beta_refined = fit_beta_binomial(filtered_bkg_nr, filtered_bkg_d, alpha_init, beta_init)
        else:
            alpha_refined = beta_refined = alpha_init = beta_init = -1

        if pid in stats:
            stats[pid][pid + bp] = {"NR" : nr , "AF" : afs, "DEPTH" : d, "ALPHA" : alpha_refined, "BETA" : beta_refined, "ALPHA_INIT" : alpha_init, "BETA_INIT" : beta_init}    
        else:
            stats[pid] = {pid + bp: {"NR" : nr , "AF" : afs, "DEPTH" : d, "ALPHA" : alpha_refined, "BETA" : beta_refined, "ALPHA_INIT" : alpha_init, "BETA_INIT" : beta_init}}

pickle.dump([stats, genes],open(outfile,"wb"))
print("Background sucessfully generated!")  