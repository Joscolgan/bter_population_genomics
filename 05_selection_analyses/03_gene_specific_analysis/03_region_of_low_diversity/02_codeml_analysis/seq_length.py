#!/usr/bin/python
from __future__ import division
from Bio import SeqIO
import sys

## Input
## python seq_length.py alignment_PAML.fasta

cmdargs = str(sys.argv)
for seq_record in SeqIO.parse(str(sys.argv[1]), "fasta"):
    output_line = '%s\t%.2f' % (seq_record.id, len(seq_record)/3.0)
    print (output_line)
