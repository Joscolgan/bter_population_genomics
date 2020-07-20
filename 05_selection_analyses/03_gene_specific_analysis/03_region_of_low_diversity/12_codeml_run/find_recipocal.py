#! /usr/bin/env python
import sys, csv

## How to use script:
## print(python script.py tblastn.output.csv blastx.output.csv output.csv)

# number below which to discard results
E_CUTOFF = 1e-3

def load_csv_to_dict(filename):
    """
    Load the CSV file into a dictionary, tying query sequences to subject
    sequences.
    """
    d = {}

    for (query_name, subject_name, score,expect) in csv.reader(open(filename)):
        # fix the names that start with a 'gi|XXXXXXX|'
        query_name = demangle_name(query_name)
        subject_name = demangle_name(subject_name)

        # convert the e-value into a floating point number
        expect = float(expect)
        if query_name not in d and expect < E_CUTOFF:
            # if we don't already have an entry for this, put it in.
            d[query_name] = subject_name

    return d

def demangle_name(name):
    """
    This functions strips off the 'gi|XXXXX|' name encoding that NCBI uses.

    Note that NCBI does this automatically for subject sequences.
    """
    if name.startswith('gi|'):
        name = name.split('|')
        name = name[2:]
        name = "|".join(name)

    return name

###

# This is the code that's run when you actually type 'find-reciprocal.py'
# at the command line; the above are function definitions, that define
# reusable blocks or chunks of code.

in_file_1 = sys.argv[1]
in_file_2 = sys.argv[2]

d1 = load_csv_to_dict(in_file_1)
d2 = load_csv_to_dict(in_file_2)

output = csv.writer(sys.stdout)

for seqname in d1:
    seqmatch1 = d1[seqname]
    seqmatch2 = d2.get(seqmatch1)

    if seqmatch2 == seqname:
        output.writerow([seqname, seqmatch1])
