#!/usr/bin/env python3
##############################################################################
##############################################################################
# Author: Joe Colgan                   Program: cnv_to_bed.py
#
# Date: 10/01/2016
#
##############################################################################
# Import modules
import sys
import re

# The purpose of the script is to take one file, reformat into BED format and output.
# The script takes two arguments from the command-line and uses the first argument
# as input and assigns the second argument as output.
# The script contains a function that does the following:
# 1) Opens the file and reads in each line.
# 2) Splits columns at tabs ('\t') and assigns the first column containing CNV type
# information (deletion/duplication) to a variable.
# 3) The second column containing information on CNV coordinates is assigned to a
# second variable.
# 4) The variables are combined and re.sub used to find and replace colons ':'
# with dashes '-'. Subsequently, the dashes are replaced with tabs ('\t').
# The formatted data is now in BED file format.
# 5) An appendable 'a' output file is opened and the contents of formatted BED file
# data are written out.

# Take input as the first command-line argument
inputfile = sys.argv[1]

# Take output as the second command-line argument
outputfile = sys.argv[2]

# Contains steps to take information from CNVnator output and convert to BED file format.
def CNV_to_BED(file_1, file_2):
    for line in open(file_1):
        first_field = line.split('\t')[0]
        second_field = line.split('\t')[1]
        coordinates = second_field + '\t' + first_field
        replaced = re.sub(':', '\t', coordinates)
        replaced_new = re.sub('-', '\t', replaced)
        outfile = open(file_2, 'a')
        outfile.write("%s\n" % replaced_new)

# Call the function
CNV_to_BED(inputfile, outputfile)
