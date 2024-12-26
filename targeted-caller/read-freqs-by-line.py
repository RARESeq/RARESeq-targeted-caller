#!/usr/bin/python3

import os
import sys
import re
import argparse

#Get inputs from command line
parser = argparse.ArgumentParser(description='Read in several freq files returning one line for each position in all the files.')
parser.set_defaults(d=':',p=0)
parser.add_argument('freqs',metavar='', help='Any number of freq files to be read in unison', nargs=argparse.REMAINDER)
parser.add_argument('-d', metavar='', type=str, help='<str> Delimiter used to separate each line. Cannot be newline character (' + str(parser.get_default('d')) + ')')
parser.add_argument('-p', metavar='', type=int, help='<int> Minimum depth. Only output positions with at least this much depth (' + str(parser.get_default('p')) + ')')
args = parser.parse_args()

#Save inputs
delimiter = args.d
depth = args.p
freqs = args.freqs

#Save reading of each file
fins = []
lines = []
for freq in freqs:
    fins.append(open(freq,'r'))
    lines.append('')

#Read in first lines
for i in range(len(lines)):
    lines[i] = fins[i].readline().strip().split('\t')
    while lines[i][0] == 'CHR': lines[i] = fins[i].readline().strip().split('\t')
#Get current chromosome as unique amongst all of them

curchr = list(set([x[0] for x in lines]))
#Make sure they all match
if len(curchr) > 1:
    print('Lines don\'t start on the same chromosome. Possibly out of order? Exiting')
    os._exit(1)
curchr = curchr[0]
#Get current position as lowest position
curpos = str(min(map(int,[x[1] for x in lines])))
#Print first line
printline = delimiter.join(['\t'.join(x) for x in lines if x[1] == curpos and int(x[2]) >= depth])
if printline != '': print(printline)

#Read remaining lines lines
while True:
    #Parse
    for i in range(len(lines)):
        #Read if it is last output positions
        if lines[i][1] == curpos and lines[i][0] == curchr:
            lines[i] = fins[i].readline().strip().split('\t')
    #Check for end of file
    while [''] in lines:
        index = lines.index([''])
        fins[index].close()
        del fins[index]
        del lines[index]
    #See of all reached the end
    if len(lines) == 0: break
    #Update chromosome if they are all in a new chromosome
    if all([x[0] != curchr for x in lines]):
        curchr = list(set([x[0] for x in lines]))
        #Make sure they all match
        if len(curchr) > 1:
            print('New chromosome out of order. Exiting')
            os._exit(1)
        curchr = curchr[0]
    #Get least position
    curpos = str(min(map(int,[x[1] for x in [y for y in lines if y[0] == curchr]])))
    printline = delimiter.join(['\t'.join(x) for x in lines if x[1] == curpos and x[0] == curchr and int(x[2]) >= depth])
    if printline != '': print(printline)