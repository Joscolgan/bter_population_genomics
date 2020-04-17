#!/usr/bin/env python3
##############################################################################
# Author: Joe Colgan                   Program: empty_BAM_filler.py
#
# Date: 07/06/2016
#
##############################################################################
## Load library system
import sys
import os

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

# Provide access to a copy of the command line arguments supplied when R session is invoked
input = sys.argv[1] # Should be file containing non-redundant list of sequences

output = sys.argv[2]

def empty_BAM_filler(inputfile, outputfile):
    #This function takes a list of files as input and for each file present in the list,\
    #evaluates if the byte size of the file isn't zero and raises an error if false
    if os.path.getsize(inputfile) == 0:
        print('%s is empty' % inputfile)
        file = open(outputfile, "w")
        file.write('%s\n' % (0))
        file.close()
    else:
        print('%s is populated' % inputfile)

empty_BAM_filler(input, output)
